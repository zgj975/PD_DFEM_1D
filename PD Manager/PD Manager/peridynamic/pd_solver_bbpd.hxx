/*==============================================================================

Copyright 2017 Dalian University of Technology .
All rights reserved

================================================================================
-- Please append file description informations here --
The Solver for 2D blank / shell
================================================================================
Date            Name                    Description of Change
2020 / 04 / 21		Zheng Guojun			Create
2020 / 04 / 21		Zheng Guojun			Update RefreshFracture
$HISTORY$
================================================================================
*/

#ifndef DLUT_SAE_PERIDYNAMIC_PD_SOLVERS_BBPD_HXX_20200421
#define DLUT_SAE_PERIDYNAMIC_PD_SOLVERS_BBPD_HXX_20200421

#include "pd_database.hxx"

namespace DLUT
{
	namespace SAE
	{
		namespace PERIDYNAMIC
		{
			//	PD Implicit analysis 
			class TSolverBBPD
			{
			public:
				TSolverBBPD() {}
				~TSolverBBPD() {}
			public:
				void				Attach(TPdModel& model)
				{
					m_pPdModel = &model;

					initializeMatParas();
				}
				/************************************************************************/
				/* 隐式求解 F=Kd 平衡方程                                               */
				/************************************************************************/
				bool				ImplicitSolve(int current_step, int total_step, int& ITERATOR_NUMS, int delta_step = 1)
				{
					TPdModel& pdModel = *m_pPdModel;
								
					ITERATOR_NUMS = 0;
					bool CONVERGENCED = FALSE;

					const set<int>& eids = pdModel.PdMeshCore().GetElementIdsByAll();
					const set<int>& nids = pdModel.PdMeshCore().GetNodeIdsByAll();
					const int nCount = (int)(nids.size());
					/************************************************************************/
					/* 节点信息初始化                                          */
					/************************************************************************/
					parallel_for_each(nids.begin(), nids.end(), [&](int nid)
						{
							TPdNode& node = pdModel.PdMeshCore().Node(nid);
							node.IteratorDisplacement().setZero();
							node.IncrementalDisplacement().setZero();

							node.Coordinate() = node.CoordinateCurrent();
						});

					parallel_for_each(eids.begin(), eids.end(), [&](int eid)
						{
							TPdElement& element = pdModel.PdMeshCore().Element(eid);
							for (int is = 0; is < IP_COUNT; ++is)
							{
								element.IP(is).Strain() = element.IP(is).StrainCurrent();
								element.IP(is).Stress() = element.IP(is).StressCurrent();
							}
						});
					/************************************************************************/
					/* P                                                           */
					/************************************************************************/
					vector<double> P;
					P.clear();
					P.resize(DIM * nCount, 0);
					for (const TLoadNodePoint& lp : pdModel.LoadNodePoints())
					{
						/************************************************************************/
						/* 静力加载的只有总力，因此将力的增量平均分布到每个增量步中             */
						/************************************************************************/
						int nid = lp.Id();

						double value = 0;
						int curid = lp.Lcid();

						if (pdModel.CurveExist(curid))
						{
							TCurve& curve = pdModel.Curve(curid);
							value = (curve.GetValueByX(current_step) - curve.GetValueByX(0)) * lp.Sf();
						}
						else
						{
							double cur_index = (double)(current_step) / (double)(total_step);
							value = lp.Sf() * cur_index;
						}

						int row = DIM * nid + (lp.Dof() - 1);
						P[row] = value;
					}

					/************************************************************************/
					/* 迭代求解                                          */
					/************************************************************************/
					while (ITERATOR_NUMS < 50 && !CONVERGENCED)
					{
						/************************************************************************/
						/* R                                          */
						/************************************************************************/
						vector<double> R;
						R.clear();
						R.resize(DIM * nCount, 0);
						for (int ni : nids)
						{
							TPdNode& node = pdModel.PdMeshCore().Node(ni);
							for (int loop_dim = 0; loop_dim < DIM; ++loop_dim)
							{
								R[ni * DIM + loop_dim] = node.InnerForce()[loop_dim];
							}
						}
						/************************************************************************/
						/* F=P-R                                          */
						/************************************************************************/
						vector<double> FORCE;
						FORCE.clear();
						FORCE.resize(DIM * nCount, 0);
						for (int loop = 0; loop < DIM * nCount; ++loop)
						{
							FORCE[loop] = P[loop] - R[loop];
						}
						/************************************************************************/
						/* 生成总刚矩阵                                                         */
						/************************************************************************/
						genGlobalStiffnessFEM();
						/************************************************************************/
						/* 施加约束条件                                                         */
						/************************************************************************/
						for (const TBoundarySpcNode& tbsn : pdModel.BoundarySpcNodes())
						{
							int nid = tbsn.Id();							
							int curSeri = DIM * nid + (tbsn.Dof() - 1);

							for (const pair<int, double>& j_v : m_GK[curSeri])
							{
								if (j_v.first != curSeri)
								{
									m_GK[j_v.first].erase(curSeri);
								}
							}

							m_GK[curSeri].clear();
							m_GK[curSeri][curSeri] = 1;

							FORCE[curSeri] = 0;
						}
						/************************************************************************/
						/* 施加强制位移                                                         */
						/************************************************************************/
						/*
						for (const TBoundaryPrescribedMotion& bpm : pdModel.BoundaryPrescribedMotions())
						{
							string type = bpm.MotionType();
							set<int>	nids;
							nids.clear();
							if (type == "RIGID")
							{
								int partId = bpm.Id();
								TPart& part = pdModel.Part(partId);
								nids = part.GetElementIds();
							}
							else if (type == "NODE")
							{
								nids.insert(bpm.Id());
							}
							for (int nid : nids)
							{
								int curid = bpm.Lcid();
								TCurve& curve = pdModel.Curve(curid);
								// 只对位移边界条件进行处理
								if (bpm.Vda() == 2)
								{
									//	按照曲线进行增量位移的计算
									double b_current = curve.GetValueByX((double)current_step / (double)total_step) * bpm.Sf();
									double b_last = curve.GetValueByX((double)(current_step - delta_step) / (double)(total_step)) * bpm.Sf();
									double b = b_current - b_last;

									if (bpm.Dof() == 3 || bpm.Dof() == 4 || bpm.Dof() == 5 || bpm.Dof() == 6)
									{
										continue;
									}

									int k = DIM * nid + (bpm.Dof() - 1);
									double Krr = m_GK[k][k];
									m_GK[k][k] = Krr * 10E10;

									FORCE[k] = Krr * b * 10E10;
								}
							}
						}
						*/
						/************************************************************************/
						/* 求解线性方程组						                                */
						/************************************************************************/
						SparseMatrix<double> sparse_matrix_GK;
						TransVecMap2SparseMatrix(m_GK, sparse_matrix_GK);

						vector<double> DISPLACEMENT;
						DISPLACEMENT.clear();
						DISPLACEMENT.resize(DIM * nCount, 0);

						umf_solver(sparse_matrix_GK, DISPLACEMENT, FORCE);
						/************************************************************************/
						/* 更新位移增量结果							                            */
						/************************************************************************/
						for (int nid : nids)
						{
							TPdNode& node = pdModel.PdMeshCore().Node(nid);
							for (int loop_dim = 0; loop_dim < DIM; ++loop_dim)
							{
								node.IteratorDisplacement()[loop_dim] = DISPLACEMENT[DIM * nid + loop_dim];
								node.IncrementalDisplacement()[loop_dim] += DISPLACEMENT[DIM * nid + loop_dim];
							}
						}
						/************************************************************************/
						/* 更新应变/应力/单元&节点内力                                          */
				 		/************************************************************************/
						updateStrainStress();
						updateInnerForce();
							
//						test();
						/************************************************************************/
						/* 迭代收敛之后，更新节点坐标和位移信息                                */
						/************************************************************************/
						CONVERGENCED = isConvergenced();
						if (CONVERGENCED)
						{											
							updateDisplacement();
							updateLocalSystem();
						
							break;
						}

						ITERATOR_NUMS++;
					}
									
					return CONVERGENCED;
				}
				/************************************************************************/
				/* 显式求解 中心差分法                                                  */
				/************************************************************************/
				void				ExplicitSolve(double time_interval)
				{
					TPdModel& pdModel = *m_pPdModel;
					calPdForces();

					for (const TPart& part : pdModel.Parts())
					{
						const TMaterial& material = pdModel.Material(part.MaterialId());
						if (material.Name() == "MAT_RIGID")
							continue;

						const TSection& section = pdModel.Section(part.SectionId());
						double Rho = material.GetMatValue("Rho");
						double thickness = section.GetSectionValue("THICKNESS");

						const set<int>& eleIds = part.GetElementIds();
						set<int> nodeIds;
						nodeIds.clear();
						for (int eid : eleIds)
						{
							const TPdElement& element = pdModel.PdMeshCore().Element(eid);
							for (int nid : element.NodeIds())
							{
								nodeIds.insert(nid);
							}
						}

						parallel_for_each(nodeIds.begin(), nodeIds.end(), [&](int nid)
						{
							TPdNode& node = pdModel.PdMeshCore().Node(nid);
							double nVolume = 0;
							set<int> adjEleIds = node.AdjElementIds();
							for (int adjEle : adjEleIds)
							{
								const TPdElement& element = pdModel.PdMeshCore().Element(adjEle);
								nVolume += 0.25 * element.Area() * element.SideLength();
							}
							node.Acceleration() = (node.OuterForce() - node.InnerForce()) / (Rho * nVolume);
							node.Velocity() = node.Velocity() + node.Acceleration() * time_interval;
							node.Displacement() = node.Displacement() + node.Velocity() * time_interval;
						});
					}
				}
				/************************************************************************/
				/* 更新断裂失效的Bond信息                                               */
				/************************************************************************/
				void				RefreshFracture()
				{
					TPdModel& pdModel = *m_pPdModel;
					for (const TPart& part : pdModel.Parts())
					{
						const TMaterial& material = pdModel.Material(part.MaterialId());
						double E = material.GetMatValue("E");
						const set<int> eids = part.GetElementIds();
	
						parallel_for_each(eids.begin(), eids.end(), [&](int ei)
						{
							TPdElement& element_i = pdModel.PdMeshCore().Element(ei);
							double s0 = element_i.CalParas().s0;
							bool updateSK = false;

							if (element_i.CalParas().b_facture)
							{
								map<int, TPdBond>& familyBonds = element_i.FamilyElementBonds();
								for (map<int, TPdBond>::iterator iter = familyBonds.begin();
									iter != familyBonds.end();)
								{
									int ej = (*iter).first;
									TPdElement& element_j = pdModel.PdMeshCore().Element(ej);
									if (element_j.CalParas().b_facture)
									{
										const TCoordinate& Xi = element_i.CoordinateInElement(0);
										const TDisplacement& Di = element_i.DisplaceInElement(0);

										const TCoordinate& Xj = element_j.CoordinateInElement(0);
										const TDisplacement& Dj = element_j.DisplaceInElement(0);

										double idist = Distance_2pt<Vector3d>(Xi, Xj);
										double nlength = Distance_2pt<Vector3d>(Xi + Di.block(0, 0, 3, 1), Xj + Dj.block(0, 0, 3, 1));
										double s = abs((nlength - idist) / idist);
										if (s > s0)
										{
											++iter;
											element_i.DeleteFamilyElement(element_j.Id());
											updateSK = true;
										}
										else
										{
											++iter;
										}									
									}
									else
									{
										++iter;
									}
								}
								/************************************************************************/
								/* 重新对单元i进行单刚的生成                                           */
								/************************************************************************/
								if (updateSK)
								{
									genSingleStiffness(element_i);
								}
							}
						});
					}
				}
			private:
				/************************************************************************/
				/* 初始化PD点所需的计算参数                                             */
				/************************************************************************/
				void				initializeMatParas()
				{
					S_IP[0] = 0.1127016653792585;
					S_IP[1] = 0.5;
					S_IP[2] = 0.8872983346207415;

					H_IP[0] = 0.277777777777778;
					H_IP[1] = 0.444444444444445;
					H_IP[2] = H_IP[0];

					TPdModel& pdModel = *m_pPdModel;
						
					for (const TPart& part : pdModel.Parts())
					{
						const TMaterial& material = pdModel.Material(part.MaterialId());
						const TSection& section = pdModel.Section(part.SectionId());
						string stype = section.Type();

						/************************************************************************/
						/*  获取材料参数                                                        */
						/************************************************************************/					
						double E, PR, Rho, G0;
						if (material.Name() == "MAT_ELASTIC")
						{
							E= material.GetMatValue("E");
							PR = material.GetMatValue("PR");
							Rho = material.GetMatValue("Rho");
							G0 = material.GetMatValue("STRESS_TENSILE");
						}
												
						double A = section.GetSectionValue("Area");
						double Iy = section.GetSectionValue("ISS");
						double Iz = section.GetSectionValue("ITT");

						D_Elastic.setZero();
						D_Elastic(0, 0) = E * A;
						D_Elastic(1, 1) = E * Iy;
						D_Elastic(2, 2) = E * Iz;
					//	/************************************************************************/
					//	/* 计算微模量                                                           */
					//	/************************************************************************/
					//	const set<int> eids = part.GetElementIds();
					//	parallel_for_each(eids.begin(), eids.end(), [&](int ei) {
					//		TPdElement& element_i = pdModel.PdMeshCore().Element(ei);

					//		const double dx = element_i.SideLength();
					//		double horizon = 0.0;
					//		if (DLUT::SAE::PERIDYNAMIC::USE_CONSTANT_HORIZON)
					//		{
					//			horizon = DLUT::SAE::PERIDYNAMIC::HORIZON;
					//		}
					//		else
					//		{
					//			horizon = element_i.SideLength() * DLUT::SAE::PERIDYNAMIC::RATIO_OF_HORIZON_MESHSIZE;
					//		}

					//		double& density = element_i.CalParas().density;
					//		double& c = element_i.CalParas().c;
					//		double& s0 = element_i.CalParas().s0;

					//		c = 9 * E / (PI * (pow(horizon, 3)));

					//		s0 = sqrt((4 * PI * G0) / (9 * E * horizon));
					//	});

					//	/************************************************************************/
					//	/* 计算所有Bond的单刚矩阵                                               */
					//	/************************************************************************/
					//	parallel_for_each(eids.begin(), eids.end(), [&](int ei)
					//	{
					//		TPdElement& element_i = pdModel.PdMeshCore().Element(ei);
					//		genSingleStiffness(element_i);
					//	});
					}
				}	
				/************************************************************************/
				/* 生成总刚矩阵		   		                                            */
				/************************************************************************/
				void				genGlobalStiffness()
				{
					TPdModel& pdModel = *m_pPdModel;

					int nCount = pdModel.PdMeshCore().NodeCount();
					int dim = 2;
					m_GK.clear();
					m_GK.resize(dim * nCount);

					for (const TPart& part : pdModel.Parts())
					{
						const set<int> eids = part.GetElementIds();
						const TSection& section = pdModel.Section(part.SectionId());
						string stype = section.Type();
							
						/************************************************************************/
						/* 组装总刚矩阵                                                         */
						/************************************************************************/
						for (int ei : eids)
						{
							const TPdElement& element_i = pdModel.PdMeshCore().Element(ei);
							//	将单元I和单元J对应的节点放入一个nids，长度为8
							const vector<int>& nids_i = element_i.NodeIds();

							const map<int, TPdBond>& familyBonds = element_i.FamilyElementBonds();
							for (const pair<int, TPdBond>& bondInfo : familyBonds)
							{
								int ej = bondInfo.first;
								const TPdElement& element_j = pdModel.PdMeshCore().Element(ej);
								const vector<int>& nids_j = element_j.NodeIds();

								vector<int> nids(nids_i);								
								for (int nj : nids_j)
								{
									nids.push_back(nj);
								}

								const SingleStiffness& SK = bondInfo.second.SK();
								int LenOfSK = 16;
							
								for (int row = 0; row < LenOfSK; ++row)
								{
									for (int col = 0; col < LenOfSK; ++col)
									{
										int Row = nids[row / dim] * dim + row % dim;
										int Col = nids[col / dim] * dim + col % dim;
										m_GK[Row][Col] += SK(row, col);
									}
								}
							}
						}
					}
				}
				/************************************************************************/
				/* 生成单元I与单元J之间的单刚矩阵                                       */
				/************************************************************************/
				void				genSingleStiffness(TPdElement& element_i)             
				{
					TPdModel& pdModel = *m_pPdModel;
					TPart& part = pdModel.Part(element_i.PartId());
					TSection& section = pdModel.Section(part.SectionId());
					double thickness = section.GetSectionValue("THICKNESS");

					map<int, TPdBond>& familyBonds = element_i.FamilyElementBonds();
					for (map<int, TPdBond>::iterator iter = familyBonds.begin();
						iter != familyBonds.end(); ++iter)
					{
						int nj = (*iter).first;
						TPdBond& bond_ij = (*iter).second;
						const TPdElement& element_j = pdModel.PdMeshCore().Element(nj);

						SingleStiffness& SK = bond_ij.SK();
						SK.resize(16, 16);
						SK.setZero();

						const double c = element_i.CalParas().c;
						double volume_scale = bond_ij.Volume() / (element_j.Area() * element_j.SideLength());

						MatrixXd k;
						k.resize(2, 2);
						k.setZero();

						k(0, 0) = 1;
						k(0, 1) = -1;
						k(1, 0) = -1;
						k(1, 1) = 1;

						MatrixXd Te = Telem(element_i, element_j);
						MatrixXd TE = TELEM(element_i, element_j);						
					}
				}
				/************************************************************************/
				/* 计算PD点之间的力密度	                                                */
				/************************************************************************/
				void				calPdForces()
				{
					TPdModel& pdModel = *m_pPdModel;

					/************************************************************************/
					/* 每一步都需要将内力置零初始化                                         */
					/************************************************************************/
					const set<int> nids = pdModel.PdMeshCore().GetNodeIdsByAll();
					parallel_for_each(nids.begin(), nids.end(), [&](int ni)
						{
							TPdNode& node = pdModel.PdMeshCore().Node(ni);
							node.InnerForce().setZero();
						});

					const set<int> eids = pdModel.PdMeshCore().GetElementIdsByAll();
					parallel_for_each(eids.begin(), eids.end(), [&](int ei)
						{
							TPdElement& element_i = pdModel.PdMeshCore().Element(ei);
							//	将单元I和单元J对应的节点放入一个nids，长度为8
							vector<int> nids_i = element_i.NodeIds();

							map<int, TPdBond>& familyBonds = element_i.FamilyElementBonds();
							for (map<int, TPdBond>::iterator iter = familyBonds.begin();
								iter != familyBonds.end(); ++iter)
							{
								int ej = iter->first;
								TPdElement& element_j = pdModel.PdMeshCore().Element(ej);
								vector<int> nids_j = element_j.NodeIds();
								const SingleStiffness& SK = iter->second.SK();

								VectorXd DIS;
								DIS.resize(16);
								DIS.setZero();
								int loop_Dis = 0;
								for (int ni : nids_i)
								{
									DIS[loop_Dis++] = pdModel.PdMeshCore().Node(ni).Displacement().x();
									DIS[loop_Dis++] = pdModel.PdMeshCore().Node(ni).Displacement().y();
								}
								for (int nj : nids_j)
								{
									DIS[loop_Dis++] = pdModel.PdMeshCore().Node(nj).Displacement().x();
									DIS[loop_Dis++] = pdModel.PdMeshCore().Node(nj).Displacement().y();
								}

								iter->second.ForceOfBond() = SK * DIS;
							}
						});

					for (int ei : eids)
					{
						TPdElement& element_i = pdModel.PdMeshCore().Element(ei);
						vector<int> nids_i = element_i.NodeIds();
						map<int, TPdBond>& familyBonds = element_i.FamilyElementBonds();
						for (const pair<int, TPdBond>& bondInfo : familyBonds)
						{
							int ej = bondInfo.first;
							TPdElement& element_j = pdModel.PdMeshCore().Element(ej);
							const VectorXd& FORCE = bondInfo.second.ForceOfBond();
							vector<int> nids_j = element_j.NodeIds();

							int loop_force = 0;
							for (int ni : nids_i)
							{
								pdModel.PdMeshCore().Node(ni).InnerForce().x() += FORCE[loop_force++];
								pdModel.PdMeshCore().Node(ni).InnerForce().y() += FORCE[loop_force++];
							}
							for (int nj : nids_j)
							{
								pdModel.PdMeshCore().Node(nj).InnerForce().x() += FORCE[loop_force++];
								pdModel.PdMeshCore().Node(nj).InnerForce().y() += FORCE[loop_force++];
							}
						}
					}
				}

			private:
				Eigen::MatrixXd N(double is, double it, double js, double jt)
				{
					Eigen::MatrixXd Res;
					Res.resize(4, 16);
					Res.setZero();
					double Ni1 = (1 - is) * (1 - it) / 4;
					double Ni2 = (1 + is) * (1 - it) / 4;
					double Ni3 = (1 + is) * (1 + it) / 4;
					double Ni4 = (1 - is) * (1 + it) / 4;
					Res(0, 0) = Res(1, 1) = Ni1;
					Res(0, 2) = Res(1, 3) = Ni2;
					Res(0, 4) = Res(1, 5) = Ni3;
					Res(0, 6) = Res(1, 7) = Ni4;

					double Nj1 = (1 - js) * (1 - jt) / 4;
					double Nj2 = (1 + js) * (1 - jt) / 4;
					double Nj3 = (1 + js) * (1 + jt) / 4;
					double Nj4 = (1 - js) * (1 + jt) / 4;
					Res(2, 8) = Res(3, 9) = Nj1;
					Res(2, 10) = Res(3, 11) = Nj2;
					Res(2, 12) = Res(3, 13) = Nj3;
					Res(2, 14) = Res(3, 15) = Nj4;

					return Res;
				}
				Eigen::MatrixXd Telem(const TPdElement& element_i, const TPdElement& element_j)
				{
					TPdModel& pdModel = *m_pPdModel;

					Eigen::MatrixXd Res;
					Res.resize(4, 4);
					Res.setZero();

					Matrix3d Local_Ei = element_i.LocalCoorSystem();
					Matrix3d Local_Ej = element_j.LocalCoorSystem();

					Res.block(0, 0, 2, 2) = Local_Ei.block(0, 0, 2, 2);
					Res.block(2, 2, 2, 2) = Local_Ej.block(0, 0, 2, 2);

					return Res;
				}
				Eigen::MatrixXd TELEM(const TPdElement& element_i, const TPdElement& element_j)
				{
					TPdModel& pdModel = *m_pPdModel;
					Eigen::MatrixXd Res;
					Res.resize(16, 16);
					Res.setZero();

					Matrix3d Local_Ei = element_i.LocalCoorSystem();
					Matrix3d Local_Ej = element_j.LocalCoorSystem();

					Res.block(0, 0, 2, 2) = Local_Ei.block(0, 0, 2, 2);
					Res.block(2, 2, 2, 2) = Local_Ei.block(0, 0, 2, 2);
					Res.block(4, 4, 2, 2) = Local_Ei.block(0, 0, 2, 2);
					Res.block(6, 6, 2, 2) = Local_Ei.block(0, 0, 2, 2);

					Res.block(8, 8, 2, 2) = Local_Ej.block(0, 0, 2, 2);
					Res.block(10, 10, 2, 2) = Local_Ej.block(0, 0, 2, 2);
					Res.block(12, 12, 2, 2) = Local_Ej.block(0, 0, 2, 2);
					Res.block(14, 14, 2, 2) = Local_Ej.block(0, 0, 2, 2);

					return Res;
				}
	
			private:
				/************************************************************************/
				/* FEM                                  */
				/************************************************************************/
				void				genSingleStiffnessFEM()
				{
					TPdModel& pdModel = *m_pPdModel;
					const set<int>& eids = pdModel.PdMeshCore().GetElementIdsByAll();
					/*parallel_*/for_each(eids.begin(), eids.end(), [&](int eid)
						{
							TPdElement& element = pdModel.PdMeshCore().Element(eid);
							SingleStiffness& SK = element.SK();
							SK.resize(12, 12);
							SK.setZero();

							double L = element.SideLength();

							SingleStiffness k;
							k.resize(12, 12);
							k.setZero();

							MatrixXd T;
							T.resize(12, 12);
							T.setZero();
							T.block(0, 0, 3, 3) = element.LocalCoorSystem();
							T.block(3, 3, 3, 3) = element.LocalCoorSystem();
							T.block(6, 6, 3, 3) = element.LocalCoorSystem();
							T.block(9, 9, 3, 3) = element.LocalCoorSystem();

							for (int is = 0; is < IP_COUNT; ++is)
							{
								Eigen::MatrixXd BL = BL_Matrix(L, S_IP[is]);
								Eigen::MatrixXd G = G_Matrix(L, S_IP[is]);

								Eigen::Matrix2d M;
								M.setZero();
								M(0, 0) = M(1, 1) = element.IP(is).StressCurrent().x();
								
								k += H_IP[is] * BL.transpose() * D_Elastic * BL * L;	// K0
								k += H_IP[is] * G.transpose() * M * G * L;				// K_sigma
							}
							k(3, 3) = k(9, 9) = D_Elastic(1, 1);
							SK = T.transpose() * k * T;	
						});
				}
				void				genGlobalStiffnessFEM()
				{
					genSingleStiffnessFEM();

					TPdModel& pdModel = *m_pPdModel;

					const set<int>& nids = pdModel.PdMeshCore().GetNodeIdsByAll();
					int nCount = (int)(nids.size());
					m_GK.clear();
					m_GK.resize(DIM * nCount);
					const set<int> eids = pdModel.PdMeshCore().GetElementIdsByAll();
					for (int ei : eids)
					{
						const TPdElement& element = pdModel.PdMeshCore().Element(ei);
						
						const vector<int>& nids = element.NodeIds();
						const SingleStiffness& SK = element.SK();
						const int LenOfSK = DIM * 2;
						for (int row = 0; row < LenOfSK; ++row)
						{
							for (int col = 0; col < LenOfSK; ++col)
							{
								int Row = nids[row / DIM] * DIM + row % DIM;
								int Col = nids[col / DIM] * DIM + col % DIM;
								m_GK[Row][Col] += SK(row, col);
							}
						}
					}
				}
				void				updateStrainStress()
				{
					TPdModel& pdModel = *m_pPdModel;
					set<int> eids = pdModel.PdMeshCore().GetElementIdsByAll();
					parallel_for_each(eids.begin(), eids.end(), [&](int eid)
						{
							TPdElement& element = pdModel.PdMeshCore().Element(eid);
							double L = element.SideLength();

							int nid1 = element.NodeId(0);
							int nid2 = element.NodeId(1);
							const TPdNode& node1 = pdModel.PdMeshCore().Node(nid1);
							const TPdNode& node2 = pdModel.PdMeshCore().Node(nid2);

							MatrixXd delta_u_global;
							delta_u_global.resize(12, 1);
							delta_u_global.block(0, 0, 6, 1) = node1.IncrementalDisplacement();
							delta_u_global.block(6, 0, 6, 1) = node2.IncrementalDisplacement();

							MatrixXd T;
							T.resize(12, 12);
							T.setZero();
							T.block(0, 0, 3, 3) = element.LocalCoorSystem();
							T.block(3, 3, 3, 3) = element.LocalCoorSystem();
							T.block(6, 6, 3, 3) = element.LocalCoorSystem();
							T.block(9, 9, 3, 3) = element.LocalCoorSystem();

							MatrixXd delta_u_local = T * delta_u_global;

							for (int is = 0; is < IP_COUNT; ++is)
							{
								Eigen::MatrixXd BL = BL_Matrix(L, S_IP[is]);
								Eigen::MatrixXd BN_star = BN_star_Matrix(L, S_IP[is], delta_u_local);

								TStrain delta_strain_local = (BL + BN_star) * delta_u_local;
								TStress delta_stress_local = D_Elastic * delta_strain_local;

								element.IP(is).StrainCurrent() = element.IP(is).Strain() + delta_strain_local;
								element.IP(is).StressCurrent() = element.IP(is).Stress() + delta_stress_local;
							}
						});
				}
				void				updateInnerForce()
				{
					TPdModel& pdModel = *m_pPdModel;
					set<int> eids = pdModel.PdMeshCore().GetElementIdsByAll();
					set<int> nids = pdModel.PdMeshCore().GetNodeIdsByAll();

					parallel_for_each(nids.begin(), nids.end(), [&](int nid)
						{
							TPdNode& node = pdModel.PdMeshCore().Node(nid);
							node.InnerForce().setZero();
						});

					for(int eid : eids)
					{
						TPdElement& element = pdModel.PdMeshCore().Element(eid);

						double L = element.SideLength();

						int nid1 = element.NodeId(0);
						int nid2 = element.NodeId(1);
						TPdNode& node1 = pdModel.PdMeshCore().Node(nid1);
						TPdNode& node2 = pdModel.PdMeshCore().Node(nid2);

						MatrixXd T;
						T.resize(12, 12);
						T.setZero();
						T.block(0, 0, 3, 3) = element.LocalCoorSystem();
						T.block(3, 3, 3, 3) = element.LocalCoorSystem();
						T.block(6, 6, 3, 3) = element.LocalCoorSystem();
						T.block(9, 9, 3, 3) = element.LocalCoorSystem();

						Eigen::MatrixXd inner_force_local;
						inner_force_local.resize(12, 1);
						inner_force_local.setZero();

						for (int is = 0; is < IP_COUNT; ++is)
						{
							Eigen::MatrixXd BL = BL_Matrix(L, S_IP[is]);

							inner_force_local += H_IP[is] * BL.transpose() * element.IP(is).StressCurrent() * L;
						}

						Eigen::MatrixXd inner_force_global = T.transpose() * inner_force_local;

						node1.InnerForce() += inner_force_global.block(0, 0, 6, 1);
						node2.InnerForce() += inner_force_global.block(6, 0, 6, 1); 
					}
				}
				void				updateLocalSystem()
				{
					TPdModel& pdModel = *m_pPdModel;
					set<int> eids = pdModel.PdMeshCore().GetElementIdsByAll();
					parallel_for_each(eids.begin(), eids.end(), [&](int eid)
						{
								TPdElement& element = pdModel.PdMeshCore().Element(eid);
									
								int nid1 = element.NodeId(0);
								int nid2 = element.NodeId(1);
								TPdNode& node1 = pdModel.PdMeshCore().Node(nid1);
								TPdNode& node2 = pdModel.PdMeshCore().Node(nid2); 
													
								Vector3d x_old = node2.Coordinate() - node1.Coordinate();
								Vector3d x_new = node2.CoordinateCurrent() - node1.CoordinateCurrent();
										
								Vector3d n = Fork_Multi<Vector3d, Vector3d>(x_old, x_new);
								Normalizer(n);
								double nx = n(0);
								double ny = n(1);
								double nz = n(2);

								double C = Calculate_COS<Vector3d>(x_old, x_new);
								double S = sqrt(1 - C * C);

								Eigen::Matrix3d R;
								R(0, 0) = nx * nx * (1 - C) + C;
								R(0, 1) = nx * ny * (1 - C) - nz * S;
								R(0, 2) = nx * nz * (1 - C) + ny * S;

								R(1, 0) = ny * nx * (1 - C) + nz * S;
								R(1, 1) = ny * ny * (1 - C) + C;
								R(1, 2) = ny * nz * (1 - C) - nx * S;

								R(2, 0) = nz * nx * (1 - C) - ny * S;
								R(2, 1) = nz * ny * (1 - C) + nx * S;
								R(2, 2) = nz * nz * (1 - C) + C;

								MatrixXd lcs_new = element.LocalCoorSystem() * R.transpose();
								element.LocalCoorSystem() = lcs_new;								
						});
				}
				void				updateDisplacement()
				{
					TPdModel& pdModel = *m_pPdModel;
					set<int> nids = pdModel.PdMeshCore().GetNodeIdsByAll();
					parallel_for_each(nids.begin(), nids.end(), [&](int nid)
						{
							TPdNode& node = pdModel.PdMeshCore().Node(nid);
							/************************************************************************/
							/* 当前增量步的位移累加到全局节点总位移                                 */
							/* 节点当前步的位移更新为上一步位移+当前增量步的位移                    */
							/************************************************************************/
							node.Displacement() += node.IncrementalDisplacement();
							node.CoordinateCurrent() = node.Coordinate() + node.IncrementalDisplacement().block(0, 0, 3, 1);
						});
				} 
				bool				isConvergenced()
				{
					bool CONVERGENCED = true;

					TPdModel& pdModel = *m_pPdModel;
					set<int> nids = pdModel.PdMeshCore().GetNodeIdsByAll();
					double MaxOfIterU = -1E20;
					double MaxOfIncrU = -1E20;
					
					for (int nid : nids)
					{
						TPdNode& node = pdModel.PdMeshCore().Node(nid);
						for (int loop_dim = 0; loop_dim < DIM; ++loop_dim)
						{
							double IterU = abs(node.IteratorDisplacement()[loop_dim]);
							double IncreU = abs(node.IncrementalDisplacement()[loop_dim]);
							if (IncreU > MaxOfIncrU)
							{
								MaxOfIncrU = IncreU;
								MaxOfIterU = IterU;
							}
						}
					}

					double ratio = MaxOfIterU / MaxOfIncrU;
					if (ratio > ALPHA_d)
					{
						CONVERGENCED = false;
//						printf("%s%5d,  %s%10.8f,  %s%10.8f,  %s%10.8f\n","FREEDOM=", MARK_OF_NODE, "Ratio=", ratio, "IterU=", MaxOfIterU, "IncrU=", MaxOfIncrU);
					}
					else
					{
						CONVERGENCED = true;
//						printf("%s%5d,  %s%10.8f,  %s%10.8f,  %s%10.8f\n\n", "FREEDOM=", MARK_OF_NODE, "Ratio=", ratio, "IterU=", MaxOfIterU, "IncrU=", MaxOfIncrU);
					}

					return CONVERGENCED;
				}
			
				Eigen::MatrixXd		G_Matrix(double L, double chi)
				{
					Eigen::MatrixXd Res;
					Res.resize(2, 12);
					Res.setZero();

					Res(0, 2) = 6 * (chi * chi - chi) / L;
					Res(0, 4) = -(3 * chi * chi - 4 * chi + 1);
					Res(0, 8) = -Res(0, 2);
					Res(0, 10) = -(3 * chi * chi - 2 * chi);

					Res(1, 1) = 6 * (chi * chi - chi) / L;
					Res(1, 5) = (3 * chi * chi - 4 * chi + 1);
					Res(1, 7) = -Res(0, 2);
					Res(1, 11) = (3 * chi * chi - 2 * chi);

					return Res;
				}
				Eigen::MatrixXd		BL_Matrix(double L, double chi)
				{
					Eigen::MatrixXd Res;
					Res.resize(3, 12);
					Res.setZero(); 
					 
					Res(0, 0) = -1.0 / L;
					Res(0, 6) = 1.0 / L;

					Res(1, 2) = (12 * chi - 6) / (L * L);
					Res(1, 4) = -(6 * chi - 4) / L;
					Res(1, 8) = Res(2, 7) = -Res(1, 2);
					Res(1, 10) = -(6 * chi - 2) / L;

					Res(2, 1) = (12 * chi - 6) / (L * L);
					Res(2, 5) = (6 * chi - 4) / L;
					Res(2, 7) = -Res(1, 2);
					Res(2, 11) = (6 * chi - 2) / L;

					return Res;
				}
				Eigen::MatrixXd		BN_star_Matrix(double L, double chi, const MatrixXd& delta_u)
				{
					Eigen::MatrixXd Res;
					Res.resize(3, 12);
					Res.setZero();
					Eigen::MatrixXd G = G_Matrix(L, chi);

					Eigen::MatrixXd delta_A;
					delta_A.resize(3, 2);
					delta_A.setZero();
					delta_A.block(0, 0, 1, 2)  = (G * delta_u).transpose(); 

					Res = 0.5 * delta_A * G;
					
					return Res;
				}
				/************************************************************************/
				/* FEM                                  */
				/************************************************************************/

				void test()
				{
					TPdModel& pdModel = *m_pPdModel;
					set<int> eids = pdModel.PdMeshCore().GetElementIdsByAll();
					set<int> nids = pdModel.PdMeshCore().GetNodeIdsByAll();

					//for (int ni : nids)
					//{
					//	const TPdNode& node = pdModel.PdMeshCore().Node(ni);
					//	cout << node.IdGlobal() << ":\t" << node.InnerForce().x() << endl;
					//}
					for (int ei : eids)
					{
						const TPdElement& element = pdModel.PdMeshCore().Element(ei);
						for (int is = 0; is < IP_COUNT; ++is)
						{
							cout << element.IdGlobal() << ":\t" << element.IP(is).StressCurrent().x() << endl;				
						}
					}					
				}
			private:
				vector< map<int, double> >	m_GK;	//	Stiffness matrix

				Eigen::Matrix3d				D_Elastic;		//	Elastic matrix
				double						S_IP[3];		//	IP index;
				double						H_IP[3];		//	IP weight coefficient;
			private:
				TPdModel*					m_pPdModel;
			};
		} // END OF PERIDYNAMIC
	} // END OF DLUT
} // END OF SAE

#endif