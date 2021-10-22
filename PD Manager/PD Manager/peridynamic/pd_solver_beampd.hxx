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

#ifndef DLUT_SAE_PERIDYNAMIC_PD_SOLVERS_BEAMBASEDPD_HXX_20210423
#define DLUT_SAE_PERIDYNAMIC_PD_SOLVERS_BEAMBASEDPD_HXX_20210423

#include "pd_database.hxx"

namespace DLUT
{
	namespace SAE
	{
		namespace PERIDYNAMIC
		{	
			//	PD Implicit analysis 
			class TSolverBeamPD
			{
			public:
				TSolverBeamPD() {}
				~TSolverBeamPD() {}
			public:
				void				Attach(TPdModel& model)
				{
					m_pPdModel = &model;

					initializeMatParasPD();
				}
				/************************************************************************/
				/* ��ʽ��� F=Kd ƽ�ⷽ��                                               */
				/************************************************************************/
				bool				ImplicitSolve(int current_step, int total_step, int& ITERATOR_NUMS, int delta_step = 1)
				{
					TPdModel& pdModel = *m_pPdModel;

 					ITERATOR_NUMS = 0;
					bool CONVERGENCED = FALSE;

					const set<int>& eids = pdModel.PdMeshCore().GetElementIdsByAll();
					const set<int>& nids = pdModel.PdMeshCore().GetNodeIdsByAll();
					const int nCount = (int)(nids.size());
					//	�ܸվ����ʼ��
					
					/************************************************************************/
					/* �ڵ���Ϣ��ʼ��												        */
					/************************************************************************/
					parallel_for_each(nids.begin(), nids.end(), [&](int nid)
						{
							TPdNode& node = pdModel.PdMeshCore().Node(nid);
							node.IteratorDisplacement().setZero();
							node.IncrementalDisplacement().setZero();
						});
					/************************************************************************/
					/* ������٤�ɽ�Ԫ�ڲ����ֵ㴦����Ϣ��ʼ��                             */
					/************************************************************************/
					parallel_for_each(eids.begin(), eids.end(), [&](int eid)
						{
							TPdElement& element = pdModel.PdMeshCore().Element(eid);
							for (int is = 0; is < IP_COUNT_1D; ++is)
							{
								element.IP(is).IteratorDisplacement().setZero();
								element.IP(is).IncrementalDisplacement().setZero();
							}
						});
					/************************************************************************/
					/* P																	*/
					/************************************************************************/
					vector<double> P;
					P.clear();
					P.resize(DOF * nCount, 0);

					for (int nid : nids)
					{
						const TPdNode& node = pdModel.PdMeshCore().Node(nid);
						const vector<TLoadNodePoint>& LPs = node.LoadNodePoint();
						for (const TLoadNodePoint& lp : LPs)
						{
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

							int row = DOF * nid + (lp.Dof() - 1);
							P[row] = value;
						}
					}
					  
					/************************************************************************/
					/* �������																*/
					/************************************************************************/
					while (ITERATOR_NUMS < MAX_INTERATOR_NUMS && !CONVERGENCED)
					{
						/************************************************************************/
						/* R																	*/
						/************************************************************************/
						vector<double> R;
						R.clear();
						R.resize(DOF * nCount, 0);
						for (int ni : nids)
						{
							TPdNode& node = pdModel.PdMeshCore().Node(ni);
							for (int loop_dim = 0; loop_dim < DOF; ++loop_dim)
							{
								R[ni * DOF + loop_dim] = node.InnerForce()[loop_dim];
							}
						}
						/************************************************************************/
						/* F=P-R                                          */
						/************************************************************************/
						vector<double> FORCE;
						FORCE.clear();
						FORCE.resize(DOF * nCount, 0);
						for (int loop = 0; loop < DOF * nCount; ++loop)
						{
							FORCE[loop] = P[loop] - R[loop];
						}
						/************************************************************************/
						/* �����ܸվ���                                                         */
						/************************************************************************/
						m_GK.clear();
						m_GK.resize(DOF * nCount);
						genGlobalStiffnessFEM();
						genGlobalStiffnessPD();
						/************************************************************************/
						/* ʩ��Լ������                                                         */
						/************************************************************************/
						for (int nid : nids)
						{
							const TPdNode& node = pdModel.PdMeshCore().Node(nid);
							const vector<TBoundarySpcNode>& BSNs = node.BoundarySpcNode();
							for (const TBoundarySpcNode& tbsn : BSNs)
							{
								int curSeri = DOF * nid + (tbsn.Dof() - 1);

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
						}
						/************************************************************************/
						/* ʩ��ǿ��λ��                                                         */
						/************************************************************************/
						for (int nid : nids)
						{
							const TPdNode& node = pdModel.PdMeshCore().Node(nid);
							const vector<TBoundaryPrescribedMotion>& BPMs = node.BoundaryPreMotion();
							for (const TBoundaryPrescribedMotion& bpm : BPMs)
							{
								int curid = bpm.Lcid();
								TCurve& curve = pdModel.Curve(curid);
								// ֻ��λ�Ʊ߽��������д���
								if (bpm.Vda() == 2)
								{
									if (ITERATOR_NUMS == 0)
									{
										//	�������߽�������λ�Ƶļ���
										double b_current = curve.GetValueByX((double)current_step / (double)total_step) * bpm.Sf();
										double b_last = curve.GetValueByX((double)(current_step - delta_step) / (double)(total_step)) * bpm.Sf();
										double b = b_current - b_last;

										int k = DOF * nid + (bpm.Dof() - 1);
										double Krr = m_GK[k][k];
										m_GK[k][k] = Krr * 10E10;

										FORCE[k] = Krr * b * 10E10;
									}
									else
									{
										int curSeri = DOF * nid + (bpm.Dof() - 1);

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
								}
							}
						}

						/************************************************************************/
						/* ������Է�����						                        		*/
						/************************************************************************/						
						vector<double> DISPLACEMENT;
						DISPLACEMENT.clear();
						DISPLACEMENT.resize(DOF * nCount, 0);

						umf_solver(m_GK, DISPLACEMENT, FORCE);
						/************************************************************************/
						/* ����λ���������							                            */
						/************************************************************************/
						for (int nid : nids)
						{
							TPdNode& node = pdModel.PdMeshCore().Node(nid);
							for (int loop_dim = 0; loop_dim < DOF; ++loop_dim)
							{
								node.IteratorDisplacement()[loop_dim] = DISPLACEMENT[DOF * nid + loop_dim];
								node.IncrementalDisplacement()[loop_dim] += DISPLACEMENT[DOF * nid + loop_dim];
							}
						}						
						/************************************************************************/
						/*  ���е�ǰ�������µ�Ӧ����Ӧ����²���                                */
						/************************************************************************/
						updateStrainStressFEM();
						updateStrainStressPD();
						parallel_for_each(nids.begin(), nids.end(), [&](int nid)
							{
								TPdNode& node = pdModel.PdMeshCore().Node(nid);
								node.InnerForce().setZero();
							});
						updateInnerForceFEM();
						updateInnerForcePD();
		
						ITERATOR_NUMS++;
						/************************************************************************/
						/* ����λ�ƽ�������ж��Ƿ��Ѿ�����                                     */
						/************************************************************************/
						CONVERGENCED = isConvergenced();
						if (CONVERGENCED)
						{
							/************************************************************************/
							/* ��������֮�󣬸��½ڵ������λ����Ϣ������������ѭ��                 */
							/************************************************************************/
							updateInfoAfterConvergenceFEM();
							updateInfoAfterConvergencePD();
						}
					}

					return CONVERGENCED;
				}
				/************************************************************************/
				/* ��ʽ��� ���Ĳ�ַ�                                                  */
				/************************************************************************/
				void				ExplicitSolve(double time_interval)
				{
					TPdModel& pdModel = *m_pPdModel;
					const set<int>& eids = pdModel.PdMeshCore().GetElementIdsByAll();
					const set<int>& nids = pdModel.PdMeshCore().GetNodeIdsByAll();
					
					//	����Ӧ�䡢Ӧ��
					updateStrainStressFEM();
					updateStrainStressPD();

					//	����ڵ�����
					parallel_for_each(nids.begin(), nids.end(), [&](int nid)
						{
							TPdNode& node = pdModel.PdMeshCore().Node(nid);
							node.InnerForce().setZero();
						});
					updateInnerForceFEM();
					updateInnerForcePD();

					//	���»��ֵ���Ϣ
					updateInfoAfterConvergenceFEM();
					updateInfoAfterConvergencePD();	
	
					//	���Ĳ�ַ����㵱ǰ��������λ��
					parallel_for_each(nids.begin(), nids.end(), [&](int nid)
						{
							TPdNode& node = pdModel.PdMeshCore().Node(nid);

							int nodePos = node.Id() * DOF;
							node.Acceleration().x() = (node.OuterForce().x() - node.InnerForce().x()) / m_GM[nodePos][nodePos];
							node.Acceleration().y() = (node.OuterForce().y() - node.InnerForce().y()) / m_GM[nodePos + 1][nodePos + 1];
							node.Acceleration().z() = (node.OuterForce().z() - node.InnerForce().z()) / m_GM[nodePos + 2][nodePos + 2];
							node.Acceleration().rx() = (node.OuterForce().rx() - node.InnerForce().rx()) / m_GM[nodePos + 3][nodePos + 3];
							node.Acceleration().ry() = (node.OuterForce().ry() - node.InnerForce().ry()) / m_GM[nodePos + 4][nodePos + 4];
							node.Acceleration().rz() = (node.OuterForce().rz() - node.InnerForce().rz()) / m_GM[nodePos + 5][nodePos + 5];

							node.Velocity() = node.Velocity() + node.Acceleration() * time_interval;
							node.IncrementalDisplacement() = node.IteratorDisplacement() = node.Velocity() * time_interval;
						});
				}
				
			private:
				/************************************************************************/
				/* Begin of PD															*/
				/************************************************************************/
				void				initializeMatParasPD()

				{
					double start = clock();

					TPdModel& pdModel = *m_pPdModel;
						
					for (const TPart& part : pdModel.Parts())
					{
						const TMaterial& material = pdModel.Material(part.MaterialId());
						const TSection& section = pdModel.Section(part.SectionId());
						string stype = section.Type();

						/************************************************************************/
						/*  ��ȡ���ϲ���                                                        */
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
						double J = Iy + Iz;

						double G = E / (2 * (1 + PR));

						D_Elastic.setZero();
						D_Elastic(0, 0) = E * A;
						D_Elastic(1, 1) = E * Iy;
						D_Elastic(2, 2) = E * Iz;
						D_Elastic(3, 3) = G * J;
						/************************************************************************/
						/* ����΢ģ��                                                           */
						/************************************************************************/
						const set<int> eids = part.GetElementIds();
						parallel_for_each(eids.begin(), eids.end(), [&](int ei) {
							TPdElement& element_i = pdModel.PdMeshCore().Element(ei);
							double alpha_i = element_i.Alpha();

							const double dx = element_i.SideLength();
							double horizon = 0.0;
							if (DLUT::SAE::PERIDYNAMIC::USE_CONSTANT_HORIZON)
							{
								horizon = DLUT::SAE::PERIDYNAMIC::HORIZON;
							}
							else
							{
								horizon = element_i.SideLength() * DLUT::SAE::PERIDYNAMIC::RATIO_OF_HORIZON_MESHSIZE;
							}

							Eigen::Matrix4d& D_FEM = element_i.CalParas().D_FEM;
							Eigen::Matrix4d& D_PD = element_i.CalParas().D_PD;
							Eigen::Matrix4d& D_H = element_i.CalParas().D_H;

							D_FEM = D_Elastic;
							double& sed_criterion = element_i.CalParas().sed_criterion;
																								
							//	ȡ��Ԫ�����ģ����ڶ������ֵ�Ϊxi
							const TIntegrationPoint& xi = element_i.IP(1);
							double modified_h2_for_dpd = 0;
							double modified_h2_for_dh = 0;

							LIST_NJ_FAMILY_ELEMENT& familyElements = element_i.FamilyElements();
							for (TPdFamilyElement& family_elem : familyElements)
							{								
								const TPdElement& element_j = pdModel.PdMeshCore().Element(family_elem.Id());
								double alpha_j = element_j.Alpha();

								double volume_scale = family_elem.VolumeIndex();
								for (int js = 0; js < IP_COUNT_1D; ++js)
								{
									const TIntegrationPoint& xj = element_j.IP(js);
									double Len = Distance_2pt(xi.Coordinate(), xj.Coordinate());

									modified_h2_for_dpd += volume_scale * Len * H_IP_1D[js] * element_j.SideLength();
									modified_h2_for_dh += volume_scale * ((alpha_i + alpha_j) / 2.0) * Len * H_IP_1D[js] * element_j.SideLength();
								}
							}

							D_PD = D_FEM * 2.0 / (A * A * modified_h2_for_dpd);
							D_H = D_FEM - D_PD * (A * A * modified_h2_for_dh) / 2.0 ;

							sed_criterion = 3.0 * G0 / (A * pow(horizon, 3));
						});						
					}

					genGlobalMass();

					double total_time = (clock() - start) / 1000;
					cout << "initializeMatParasPD():\t\t" << total_time << endl;
				}	
				void				genSingleStiffnessPD()             
				{
					TPdModel& pdModel = *m_pPdModel;
					set<int> eids = pdModel.PdMeshCore().GetElementIdsByAll();
					parallel_for_each(eids.begin(), eids.end(), [&](int ei)
						{
							TPdElement& element_i = pdModel.PdMeshCore().Element(ei);
							//	�ڴ�FEM������Ҫ����PD�ĸն������
							if (element_i.AnalysisElementType() != FEM_ELEMENT)
							{
								LIST_NJ_FAMILY_ELEMENT& familyElements = element_i.FamilyElements();
								for (TPdFamilyElement& family_elem : familyElements)
								{									
									genSingleStiffnessPD_Element(ei, family_elem);
								}
							}
						});	
				}
				void				genSingleStiffnessPD_Element(int ei, TPdFamilyElement& family_elem)
				{
					TPdModel& pdModel = *m_pPdModel;
					TPdElement& element_i = pdModel.PdMeshCore().Element(ei);
					Eigen::Matrix4d& D_PD = element_i.CalParas().D_PD;

					double alpha_i = element_i.Alpha();

					const TPdElement& element_j = pdModel.PdMeshCore().Element(family_elem.Id());
					double alpha_j = element_j.Alpha();

					double alpha = (alpha_i + alpha_j) / 2.0;

					MATRIX_SINGLE_STIFFNESS_PD& SK = family_elem.SK();
					SK.setZero();

					double volume_scale = family_elem.VolumeIndex();

					//	����ÿһ��bond�ĵ���
					const vector<TPdBond>& bonds = family_elem.Bonds();
					for (const TPdBond& bond : bonds)
					{
						//	ֻ���δʧЧ��bond���в���
						if (bond.IsValid())
						{
							double L = bond.BondLength();
							const TBondPoint& Xi = bond.Xi();
							const TBondPoint& Xj = bond.Xj();

							Eigen::Matrix<double, 12, 12> k;
							k.setZero();
							for (int is = 0; is < IP_COUNT_1D; ++is)
							{
								Eigen::MatrixXd BL = BL_Matrix(L, S_IP_1D[is]);
								Eigen::MatrixXd G = G_Matrix(L, S_IP_1D[is]);

								Eigen::Matrix2d M;
								M.setZero();
								M(0, 0) = M(1, 1) = bond.IP(is).StressCurrent().x();

								k += H_IP_1D[is] * BL.transpose() * element_i.CalParas().D_PD * BL * L;	// K0
								k += H_IP_1D[is] * G.transpose() * M * G * L;							// K_sigma
							}

							Eigen::Matrix<double, 12, 24> TbTetNijTE = Tb_Tet_Nij_TE(element_i, element_j, bond);
							SK += 0.5 * alpha * volume_scale * element_i.SideLength() * element_i.Area() * H_IP_1D[Xj.Index()] * element_j.SideLength() * element_j.Area() * TbTetNijTE.transpose() * k * TbTetNijTE;
						}
					}
				}
				void				genGlobalStiffnessPD()
				{
					genSingleStiffnessPD();

					TPdModel& pdModel = *m_pPdModel;
					const set<int>& nids = pdModel.PdMeshCore().GetNodeIdsByAll();
					const set<int>& eids = pdModel.PdMeshCore().GetElementIdsByAll();
					/************************************************************************/
					/* ��װ�ܸվ���                                                         */
					/************************************************************************/
					for (int ei : eids)
					{
						const TPdElement& element_i = pdModel.PdMeshCore().Element(ei);
						if (element_i.AnalysisElementType() == PD_ELEMENT ||
							element_i.AnalysisElementType() == MORPHING_ELEMENT)
						{
							//	����ԪI�͵�ԪJ��Ӧ�Ľڵ����һ��nids������Ϊ8
							const vector<int>& nids_i = element_i.NodeIds();

							const LIST_NJ_FAMILY_ELEMENT& familyElements = element_i.FamilyElements();
							for (const TPdFamilyElement& family_elem : familyElements)
							{
								const TPdElement& element_j = pdModel.PdMeshCore().Element(family_elem.Id());
								const vector<int>& nids_j = element_j.NodeIds();

								vector<int> nids(nids_i);
								for (int nj : nids_j)
								{
									nids.push_back(nj);
								}

								const MATRIX_SINGLE_STIFFNESS_PD& SK = family_elem.SK();
								int LenOfSK = 24;

								for (int row = 0; row < LenOfSK; ++row)
								{
									for (int col = 0; col < LenOfSK; ++col)
									{
										int Row = nids[row / DOF] * DOF + row % DOF;
										int Col = nids[col / DOF] * DOF + col % DOF;
										m_GK[Row][Col] += SK(row, col);
									}
								}
							}
						}
					}					
				}
				void				updateStrainStressPD()
				{
					TPdModel& pdModel = *m_pPdModel;
					set<int> eids = pdModel.PdMeshCore().GetElementIdsByAll();
					//	���·�����٤�ɽ�Ԫ���ֵ㴦��λ��������Ϣ
					parallel_for_each(eids.begin(), eids.end(), [&](int eid)
						{
							TPdElement& element = pdModel.PdMeshCore().Element(eid);
							int nid1 = element.NodeId(0);
							int nid2 = element.NodeId(1);
							const TPdNode& node1 = pdModel.PdMeshCore().Node(nid1);
							const TPdNode& node2 = pdModel.PdMeshCore().Node(nid2);
							MatrixXd delta_u_global;
							delta_u_global.resize(12, 1);
							delta_u_global.block(0, 0, 6, 1) = node1.IteratorDisplacement();
							delta_u_global.block(6, 0, 6, 1) = node2.IteratorDisplacement();
							
							MatrixXd TE;
							TE.resize(12, 12);
							TE.setZero();
							TE.block(0, 0, 3, 3) = element.LocalCoorSystem();
							TE.block(3, 3, 3, 3) = element.LocalCoorSystem();
							TE.block(6, 6, 3, 3) = element.LocalCoorSystem();
							TE.block(9, 9, 3, 3) = element.LocalCoorSystem();

							MatrixXd TN;
							TN.resize(6, 6);
							TN.setZero();
							TN.block(0, 0, 3, 3) = element.LocalCoorSystem();
							TN.block(3, 3, 3, 3) = element.LocalCoorSystem();

							MatrixXd delta_u_local = TE * delta_u_global;

							double L = element.SideLength();
							//	���µ�Ԫ���ֵ㴦��λ��������ȫ������ϵ��
							for (int is = 0; is < IP_COUNT_1D; ++is)
							{
								Eigen::MatrixXd N = N_SF_BEAM(L, S_IP_1D[is]);
								//	ͨ��TN���ֲ�����ϵ�µĽڵ����������������ȫ������ϵ��
								element.IP(is).IteratorDisplacement() = TN.transpose() * N * delta_u_local;
								element.IP(is).IncrementalDisplacement() += element.IP(is).IteratorDisplacement();
							}
						});

					//	����bond���ֵ㴦��Ӧ���Ӧ��
					parallel_for_each(eids.begin(), eids.end(), [&](int ei)
					{
							TPdElement& element_i = pdModel.PdMeshCore().Element(ei);
							//	�Դ�FEM����������bond�������Ϣ
							if (element_i.AnalysisElementType() != FEM_ELEMENT)
							{
								double alpha_i = element_i.Alpha();
								const Eigen::Matrix4d& D_PD = element_i.CalParas().D_PD;
								LIST_NJ_FAMILY_ELEMENT& familyElements = element_i.FamilyElements();
								for (TPdFamilyElement& family_elem : familyElements)
								{									
									const TPdElement& element_j = pdModel.PdMeshCore().Element(family_elem.Id());
									double alpha_j = element_j.Alpha();
									double alpha = (alpha_i + alpha_j) / 2.0;
									double volume_scale = family_elem.VolumeIndex();

									vector<TPdBond>& bonds = family_elem.Bonds();
									for (TPdBond& bond : bonds)
									{
										const TBondPoint& Xi = bond.Xi();
										const TBondPoint& Xj = bond.Xj();

										double L = bond.BondLength();

										MatrixXd delta_u_global;
										delta_u_global.resize(12, 1);
										delta_u_global.block(0, 0, 6, 1) = Xi.IteratorDisplacement();
										delta_u_global.block(6, 0, 6, 1) = Xj.IteratorDisplacement();

										MatrixXd T;
										T.resize(12, 12);
										T.setZero();
										T.block(0, 0, 3, 3) = bond.LocalCoorSystem();
										T.block(3, 3, 3, 3) = bond.LocalCoorSystem();
										T.block(6, 6, 3, 3) = bond.LocalCoorSystem();
										T.block(9, 9, 3, 3) = bond.LocalCoorSystem();

										MatrixXd delta_u_local = T * delta_u_global;
										for (int js = 0; js < IP_COUNT_1D; ++js)
										{
											Eigen::MatrixXd BL = BL_Matrix(L, S_IP_1D[js]);
											Eigen::MatrixXd BN_star = BN_star_Matrix(L, S_IP_1D[js], delta_u_local);

											TStrain delta_strain_local = (BL + BN_star) * delta_u_local;
											TStress delta_stress_local = D_PD * delta_strain_local;
										
											bond.MicroPotentialDensity() += ((bond.IP(js).StressCurrent() + delta_stress_local * 0.5).transpose() * delta_strain_local)(0, 0) * H_IP_1D[js];

											bond.IP(js).StrainCurrent() += delta_strain_local;
											bond.IP(js).StressCurrent() += delta_stress_local;
										}
									}
								}
							}
					});
				}
				void				updateInnerForcePD()
				{
					TPdModel& pdModel = *m_pPdModel;
					set<int> eids = pdModel.PdMeshCore().GetElementIdsByAll();
					set<int> nids = pdModel.PdMeshCore().GetNodeIdsByAll();
								
					//	����bond�Ľṹ����
					for(int ei : eids)
					{
						TPdElement& element_i = pdModel.PdMeshCore().Element(ei);
						if (element_i.AnalysisElementType() == FEM_ELEMENT)
						{
							continue;
						}

						double alpha_i = element_i.Alpha();

						TPdNode& node1 = pdModel.PdMeshCore().Node(element_i.NodeId(0));
						TPdNode& node2 = pdModel.PdMeshCore().Node(element_i.NodeId(1));

						LIST_NJ_FAMILY_ELEMENT& familyElements = element_i.FamilyElements();
						for (TPdFamilyElement& family_elem : familyElements)
						{
							const TPdElement& element_j = pdModel.PdMeshCore().Element(family_elem.Id());
							double alpha_j = element_j.Alpha();

							double alpha = (alpha_i + alpha_j) / 2.0;

							TPdNode& node3 = pdModel.PdMeshCore().Node(element_j.NodeId(0));
							TPdNode& node4 = pdModel.PdMeshCore().Node(element_j.NodeId(1));

							double volume_scale = family_elem.VolumeIndex();

							vector<TPdBond>& bonds = family_elem.Bonds();
							for (TPdBond& bond : bonds)
							{
								const TBondPoint& Xi = bond.Xi();
								const TBondPoint& Xj = bond.Xj();

								double L = bond.BondLength();

								Eigen::MatrixXd r;
								r.resize(12, 1);
								r.setZero();

								for (int is = 0; is < IP_COUNT_1D; ++is)
								{
									Eigen::MatrixXd BL = BL_Matrix(L, S_IP_1D[is]);
									r += H_IP_1D[is] * BL.transpose() * bond.IP(is).StressCurrent() * L;
								}
								Eigen::Matrix<double, 12, 24> TbTetNijTE = Tb_Tet_Nij_TE(element_i, element_j, bond);
								Eigen::MatrixXd inner_force = 0.5 * alpha * volume_scale * element_i.SideLength() * element_i.Area() * H_IP_1D[Xj.Index()] * element_j.SideLength() * element_j.Area() * TbTetNijTE.transpose() * r;

								//	���ṹ�����ַ���������Ԫ���ĸ��ڵ���ȥ
								node1.InnerForce() += inner_force.block(0, 0, 6, 1);
								node2.InnerForce() += inner_force.block(6, 0, 6, 1);
								node3.InnerForce() += inner_force.block(12, 0, 6, 1);
								node4.InnerForce() += inner_force.block(18, 0, 6, 1);
							}
						}
					}
				}
				
				void				updateInfoAfterConvergencePD()
				{
					TPdModel& pdModel = *m_pPdModel;
					set<int> eids = pdModel.PdMeshCore().GetElementIdsByAll();

					parallel_for_each(eids.begin(), eids.end(), [&](int ei)
						{
							TPdElement& element_i = pdModel.PdMeshCore().Element(ei);
			
							if (element_i.AnalysisElementType() != FEM_ELEMENT)
							{
								const Eigen::Matrix4d& D_PD = element_i.CalParas().D_PD;
								TPdNode& node1 = pdModel.PdMeshCore().Node(element_i.NodeId(0));
								TPdNode& node2 = pdModel.PdMeshCore().Node(element_i.NodeId(1));

								LIST_NJ_FAMILY_ELEMENT& familyElements = element_i.FamilyElements();
								for (TPdFamilyElement& family_elem : familyElements)
								{
									const TPdElement& element_j = pdModel.PdMeshCore().Element(family_elem.Id());

									vector<TPdBond>& bonds = family_elem.Bonds();
									for (TPdBond& bond : bonds)
									{
										const TBondPoint& Xi = bond.Xi();
										const TBondPoint& Xj = bond.Xj();
										//	����bond����ϵ
										{
											Vector3d x_old = Xj.Coordinate() - Xi.Coordinate();
											Normalizer(x_old);

											TCoordinate n2_coord_new = Xj.Coordinate() + Xj.IncrementalDisplacement().block(0, 0, 3, 1);
											TCoordinate n1_coord_new = Xi.Coordinate() + Xi.IncrementalDisplacement().block(0, 0, 3, 1);
											Vector3d x_new = n2_coord_new - n1_coord_new;
											Normalizer(x_new);

											Vector3d z_new = Fork_Multi<Vector3d, Vector3d>(x_old, x_new);
											Normalizer(z_new);
											double angle_z = acos(Calculate_COS<Vector3d>(x_old, x_new));
											Eigen::Matrix3d R = RotationMatrix(z_new, angle_z);

											Vector3d angle_in_global_sys = (Xi.IncrementalDisplacement().block(3, 0, 3, 1) + Xj.IncrementalDisplacement().block(3, 0, 3, 1)) / 2.0;
											Vector3d angle_in_local_sys = bond.LocalCoorSystem() * angle_in_global_sys;

											double angle_torison = (node1.IncrementalDisplacement()(3) + node2.IncrementalDisplacement()(3)) / 2.0;
											Eigen::Matrix3d R_Torison = RotationMatrix(x_old, angle_in_local_sys(0));

											MatrixXd lcs_new = bond.LocalCoorSystem() * R_Torison.transpose() * R.transpose();
											bond.LocalCoorSystem() = lcs_new;
										}

										//	����bond���ֵ㴦΢���ܡ�Ӧ�䡢Ӧ������Ϣ
										{
											double L = bond.BondLength();

											MatrixXd delta_u_global;
											delta_u_global.resize(12, 1);
											delta_u_global.block(0, 0, 6, 1) = Xi.IncrementalDisplacement();
											delta_u_global.block(6, 0, 6, 1) = Xj.IncrementalDisplacement();

											MatrixXd T;
											T.resize(12, 12);
											T.setZero();
											T.block(0, 0, 3, 3) = bond.LocalCoorSystem();
											T.block(3, 3, 3, 3) = bond.LocalCoorSystem();
											T.block(6, 6, 3, 3) = bond.LocalCoorSystem();
											T.block(9, 9, 3, 3) = bond.LocalCoorSystem();

											MatrixXd delta_u_local = T * delta_u_global;
											for (int js = 0; js < IP_COUNT_1D; ++js)
											{
												Eigen::MatrixXd BL = BL_Matrix(L, S_IP_1D[js]);
												Eigen::MatrixXd BN_star = BN_star_Matrix(L, S_IP_1D[js], delta_u_local);

												TStrain delta_strain_local = (BL + BN_star) * delta_u_local;
												TStress delta_stress_local = D_PD * delta_strain_local;

												bond.MicroPotentialDensity() += ((bond.IP(js).StressLaststep() + delta_stress_local * 0.5).transpose() * delta_strain_local)(0, 0) * H_IP_1D[js];

												bond.IP(js).StrainLaststep() = bond.IP(js).StrainCurrent();
												bond.IP(js).StressLaststep() = bond.IP(js).StressCurrent();
											}
										}
									}
								}
							}
						});

					//	����bond�����˵㴦��������Ϣ����FEM��Ԫ�Ļ��ֵ㴦��������Ϣ����
					parallel_for_each(eids.begin(), eids.end(), [&](int eid)
						{
							TPdElement& element = pdModel.PdMeshCore().Element(eid);
							int nid1 = element.NodeId(0);
							int nid2 = element.NodeId(1);
							TPdNode& node1 = pdModel.PdMeshCore().Node(nid1);
							TPdNode& node2 = pdModel.PdMeshCore().Node(nid2);

							Eigen::VectorXd coord;
							coord.resize(6);
							coord(0) = node1.Coordinate().x();
							coord(1) = node1.Coordinate().y();
							coord(2) = node1.Coordinate().z();

							coord(3) = node2.Coordinate().x();
							coord(4) = node2.Coordinate().y();
							coord(5) = node2.Coordinate().z();

							for (int is = 0; is < IP_COUNT_1D; ++is)
							{
								//	������ֵ㴦������ʱ�ٶ����κ������Ȼ����ֱ����ʽ
								Eigen::MatrixXd N = N_SF_ROD(S_IP_1D[is]);
								element.IP(is).Coordinate() = N * coord;
							}
						});
				}
			private: 
				void				calPdForces()
				{
					TPdModel& pdModel = *m_pPdModel;
					/************************************************************************/
					/* ÿһ������Ҫ�����������ʼ��                                         */
					/************************************************************************/
					const set<int> nids = pdModel.PdMeshCore().GetNodeIdsByAll();
					parallel_for_each(nids.begin(), nids.end(), [&](int ni)
						{
							TPdNode& node = pdModel.PdMeshCore().Node(ni);
							node.InnerForce().setZero();
						});

					/************************************************************************/
					/* ���㵥ԪI�뵥ԪJ֮��������� 24*1 ����                              */
					/************************************************************************/
					const set<int> eids = pdModel.PdMeshCore().GetElementIdsByAll();
					parallel_for_each(eids.begin(), eids.end(), [&](int ei)
						{
							TPdElement& element_i = pdModel.PdMeshCore().Element(ei);

							LIST_NJ_FAMILY_ELEMENT& familyElements = element_i.FamilyElements();
							for (TPdFamilyElement& family_elem : familyElements)
							{
								TPdElement& element_j = pdModel.PdMeshCore().Element(family_elem.Id());
						
								//	����ԪI�͵�ԪJ��Ӧ�Ľڵ����һ��nids������Ϊ4
								vector<int> nids = element_i.NodeIds();
								for (int nj : element_j.NodeIds())
								{
									nids.push_back(nj);
								}
								const MATRIX_SINGLE_STIFFNESS_PD& SK = family_elem.SK();

								VectorXd DIS;
								DIS.resize(24);
								DIS.setZero();
								int loop_Dis = 0;
								for (int nid : nids)
								{
									DIS.block(6 * loop_Dis++, 0, 6, 1) = pdModel.PdMeshCore().Node(nid).TotalDisplacement();
								}

								family_elem.ForceOfBond() = SK * DIS;
							}
						});

					/************************************************************************/
					/* �����е��������ۼӵ���Ԫ�ڵ���ȥ                                     */
					/************************************************************************/
					for (int ei : eids)
					{
						TPdElement& element_i = pdModel.PdMeshCore().Element(ei);						
						LIST_NJ_FAMILY_ELEMENT& familyElements = element_i.FamilyElements();
						for (TPdFamilyElement& family_elem : familyElements)
						{
							TPdElement& element_j = pdModel.PdMeshCore().Element(family_elem.Id());
							//	����ԪI�͵�ԪJ��Ӧ�Ľڵ����һ��nids������Ϊ4
							vector<int> nids = element_i.NodeIds();
							for (int nj : element_j.NodeIds())
							{
								nids.push_back(nj);
							}
							const MATRIX_PD_FORCE_PD& FORCE = family_elem.ForceOfBond();
							
							int loop_force = 0;
							for (int ni : nids)
							{
								pdModel.PdMeshCore().Node(ni).InnerForce() += FORCE.block(6 * loop_force++, 0, 6, 1);
							}
						}
					}
				}
			private:
				Eigen::Matrix<double, 12, 24> Tb_Tet_Nij_TE(const TPdElement& element_i, const TPdElement& element_j, const TPdBond& bond)
				{
					Eigen::Matrix<double, 12, 24> Res;
					Res.setZero();

					const TBondPoint& Xi = bond.Xi();
					const TBondPoint& Xj = bond.Xj();

					Eigen::Matrix<double, 6, 6> lamda;
					lamda.setZero();
					lamda.block(0, 0, 3, 3) = bond.LocalCoorSystem();
					lamda.block(3, 3, 3, 3) = bond.LocalCoorSystem();

					Eigen::Matrix<double, 6, 6> Tei;
					Tei.setZero();
					Tei.block(0, 0, 3, 3) = element_i.LocalCoorSystem();
					Tei.block(3, 3, 3, 3) = element_i.LocalCoorSystem();

					Matrix<double, 6, 12> Ni = N_SF_BEAM(element_i.SideLength(), S_IP_1D[Xi.Index()]);

					Eigen::Matrix<double, 12, 12> TEi;
					TEi.setZero();
					TEi.block(0, 0, 6, 6) = Tei;
					TEi.block(6, 6, 6, 6) = Tei;

					Res.block(0, 0, 6, 12) = lamda * Tei.transpose() * Ni * TEi;

					Eigen::Matrix<double, 6, 6> Tej;
					Tej.setZero();
					Tej.block(0, 0, 3, 3) = element_j.LocalCoorSystem();
					Tej.block(3, 3, 3, 3) = element_j.LocalCoorSystem();

					Matrix<double, 6, 12> Nj = N_SF_BEAM(element_j.SideLength(), S_IP_1D[Xj.Index()]);

					Eigen::Matrix<double, 12, 12> TEj;
					TEj.setZero();
					TEj.block(0, 0, 6, 6) = Tej;
					TEj.block(6, 6, 6, 6) = Tej;

					Res.block(6, 12, 6, 12) = lamda * Tej.transpose() * Nj * TEj;

					return Res;
				}
				/************************************************************************/
				/* ���¶���ʧЧ��Bond��Ϣ                                               */
				/************************************************************************/
			public:
				void				RefreshFracture()
				{
					TPdModel& pdModel = *m_pPdModel;
					set<int> eids = pdModel.PdMeshCore().GetElementIdsByAll();

					parallel_for_each(eids.begin(), eids.end(), [&](int ei)
						{
							TPdElement& element_i = pdModel.PdMeshCore().Element(ei);
							// ֻ�ܶԴ�PD�����ж��Ѳ���
							if (element_i.AnalysisElementType() == PD_ELEMENT)
							{
								double sed_criterion = element_i.CalParas().sed_criterion;

								LIST_NJ_FAMILY_ELEMENT& familyElements = element_i.FamilyElements();

								for (TPdFamilyElement& family_elem : familyElements)
								{
									const TPdElement& element_j = pdModel.PdMeshCore().Element(family_elem.Id());

									if (element_i.Id() == element_j.Id())
										continue;

									vector<TPdBond>& bonds = family_elem.Bonds();
									int nCountBondFractured = 0;
									for (TPdBond& bond : bonds)
									{
										if (bond.MicroPotentialDensity() > sed_criterion)
										{
											nCountBondFractured++;
										}
									}
									if (nCountBondFractured == (int)(bonds.size()))
									{
										for (TPdBond& bond :bonds)
										{
											bond.Failured();
										}
									}
								}
							}
						});
				}
				/************************************************************************/
				/* End of PD                                  							*/
				/************************************************************************/
			private:
				/************************************************************************/
				/* Begin of FEM															*/
				/************************************************************************/

				void				genSingleStiffnessFEM()
				{
					TPdModel& pdModel = *m_pPdModel;
					const set<int>& eids = pdModel.PdMeshCore().GetElementIdsByAll();
					parallel_for_each(eids.begin(), eids.end(), [&](int eid)
						{
							TPdElement& element = pdModel.PdMeshCore().Element(eid);
							//	�ڴ�PD������Ҫ��������Ԫ���ռ���
							if (element.AnalysisElementType() != PD_ELEMENT)
							{
								MATRIX_SINGLE_STIFFNESS_FEM& SK = element.SK();
								SK.setZero();

								const Eigen::Matrix4d& D_H = element.CalParas().D_H;

								double L = element.SideLength();

								MATRIX_SINGLE_STIFFNESS_FEM k;
								k.setZero();

								MatrixXd T;
								T.resize(12, 12);
								T.setZero();
								T.block(0, 0, 3, 3) = element.LocalCoorSystem();
								T.block(3, 3, 3, 3) = element.LocalCoorSystem();
								T.block(6, 6, 3, 3) = element.LocalCoorSystem();
								T.block(9, 9, 3, 3) = element.LocalCoorSystem();

								for (int is = 0; is < IP_COUNT_1D; ++is)
								{
									Eigen::MatrixXd BL = BL_Matrix(L, S_IP_1D[is]);
									Eigen::MatrixXd G = G_Matrix(L, S_IP_1D[is]);

									Eigen::Matrix2d M;
									M.setZero();
									M(0, 0) = M(1, 1) = element.IP(is).StressCurrent().x();

									k += H_IP_1D[is] * BL.transpose() * D_H * BL * L;			// K0
									k += H_IP_1D[is] * G.transpose() * M * G * L;				// K_sigma
								}
								SK = T.transpose() * k * T;
							}
						});
				}
				void				genGlobalStiffnessFEM()
				{
					genSingleStiffnessFEM();

					TPdModel& pdModel = *m_pPdModel;

					const set<int>& nids = pdModel.PdMeshCore().GetNodeIdsByAll();
					
					const set<int> eids = pdModel.PdMeshCore().GetElementIdsByAll();
					for (int ei : eids)
					{
						const TPdElement& element = pdModel.PdMeshCore().Element(ei);
						if (element.AnalysisElementType() == FEM_ELEMENT ||
							element.AnalysisElementType() == MORPHING_ELEMENT)
						{
							const vector<int>& nids = element.NodeIds();
							const MATRIX_SINGLE_STIFFNESS_FEM& SK = element.SK();
							const int LenOfSK = DOF * 2;
							for (int row = 0; row < LenOfSK; ++row)
							{
								for (int col = 0; col < LenOfSK; ++col)
								{
									int Row = nids[row / DOF] * DOF + row % DOF;
									int Col = nids[col / DOF] * DOF + col % DOF;
									m_GK[Row][Col] += SK(row, col);
								}
							}
						}
					}
				}
				void				updateStrainStressFEM()
				{
					TPdModel& pdModel = *m_pPdModel;
					set<int> eids = pdModel.PdMeshCore().GetElementIdsByAll();
					parallel_for_each(eids.begin(), eids.end(), [&](int eid)
						{
							TPdElement&  element = pdModel.PdMeshCore().Element(eid);
							//	�ڴ�PD������Ҫ����PD�ĸն������
							if (element.AnalysisElementType() != PD_ELEMENT)
							{
								double L = element.SideLength();
								const Eigen::Matrix4d& D_H = element.CalParas().D_H;
								const Eigen::Matrix4d& D = element.CalParas().D_FEM;

								int nid1 = element.NodeId(0);
								int nid2 = element.NodeId(1);
								const TPdNode& node1 = pdModel.PdMeshCore().Node(nid1);
								const TPdNode& node2 = pdModel.PdMeshCore().Node(nid2);

								MatrixXd delta_u_global;
								delta_u_global.resize(12, 1);
								delta_u_global.block(0, 0, 6, 1) = node1.IteratorDisplacement();
								delta_u_global.block(6, 0, 6, 1) = node2.IteratorDisplacement();

								MatrixXd T;
								T.resize(12, 12);
								T.setZero();
								T.block(0, 0, 3, 3) = element.LocalCoorSystem();
								T.block(3, 3, 3, 3) = element.LocalCoorSystem();
								T.block(6, 6, 3, 3) = element.LocalCoorSystem();
								T.block(9, 9, 3, 3) = element.LocalCoorSystem();

								MatrixXd delta_u_local = T * delta_u_global;

								for (int is = 0; is < IP_COUNT_1D; ++is)
								{
									Eigen::MatrixXd BL = BL_Matrix(L, S_IP_1D[is]);
									Eigen::MatrixXd BN_star = BN_star_Matrix(L, S_IP_1D[is], delta_u_local);

									TStrain delta_strain_local = (BL + BN_star) * delta_u_local;
									TStress delta_stress_local = D_H * delta_strain_local;
																		
									element.IP(is).StrainCurrent() += delta_strain_local;
									element.IP(is).StressCurrent() += delta_stress_local;
								}
							}							
						});
				}
				void				updateInnerForceFEM()
				{
					TPdModel& pdModel = *m_pPdModel;
					set<int> eids = pdModel.PdMeshCore().GetElementIdsByAll();
					set<int> nids = pdModel.PdMeshCore().GetNodeIdsByAll();

					for(int eid : eids)
					{
						TPdElement& element = pdModel.PdMeshCore().Element(eid);
						//	�ڴ�PD������Ҫ����PD�ĸն������
						if (element.AnalysisElementType() != PD_ELEMENT)
						{
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

							for (int is = 0; is < IP_COUNT_1D; ++is)
							{
								Eigen::MatrixXd BL = BL_Matrix(L, S_IP_1D[is]);

								inner_force_local += H_IP_1D[is] * BL.transpose() * element.IP(is).StressCurrent() * L;
							}

							Eigen::MatrixXd inner_force_global = T.transpose() * inner_force_local;

							node1.InnerForce() += inner_force_global.block(0, 0, 6, 1);
							node2.InnerForce() += inner_force_global.block(6, 0, 6, 1);
						}
					}
				}
				void				updateInfoAfterConvergenceFEM()
				{
					TPdModel& pdModel = *m_pPdModel;
					set<int> eids = pdModel.PdMeshCore().GetElementIdsByAll();
					set<int> nids = pdModel.PdMeshCore().GetNodeIdsByAll();
			
					//	����FEM��Ԫ��Ϣ
					parallel_for_each(eids.begin(), eids.end(), [&](int eid)
						{
							TPdElement& element = pdModel.PdMeshCore().Element(eid);
							int nid1 = element.NodeId(0);
							int nid2 = element.NodeId(1);
							TPdNode& node1 = pdModel.PdMeshCore().Node(nid1);
							TPdNode& node2 = pdModel.PdMeshCore().Node(nid2);

							//	���µ�Ԫ����ϵ
							{
								Vector3d x_old = node2.Coordinate() - node1.Coordinate();
								Normalizer(x_old);
								TCoordinate n2_coord_new = node2.Coordinate() + node2.IncrementalDisplacement().block(0, 0, 3, 1);
								TCoordinate n1_coord_new = node1.Coordinate() + node1.IncrementalDisplacement().block(0, 0, 3, 1);

								Vector3d x_new = n2_coord_new - n1_coord_new;
								Normalizer(x_new);

								double angle_z = acos(Calculate_COS<Vector3d>(x_old, x_new));
								//	����Ƕ�Լ����0�������κδ����������ֵ�ǰ����ϵ
								if (angle_z > ERR_VALUE)
								{
									Vector3d z_new = Fork_Multi<Vector3d, Vector3d>(x_old, x_new);
									Normalizer(z_new);
									Eigen::Matrix3d R_RotationOfZ = RotationMatrix(z_new, angle_z);

									Vector3d angle_in_global_sys = (node1.IncrementalDisplacement().block(3, 0, 3, 1) + node2.IncrementalDisplacement().block(3, 0, 3, 1)) / 2.0;
									Vector3d angle_in_local_sys = element.LocalCoorSystem() * angle_in_global_sys;
									Eigen::Matrix3d R_TorisonOfX = RotationMatrix(x_old, angle_in_local_sys(0));

									MatrixXd lcs_new = element.LocalCoorSystem() * R_TorisonOfX.transpose() * R_RotationOfZ.transpose();
									element.LocalCoorSystem() = lcs_new;
								}
							}							

							//	���µ�Ԫ���ֵ㴦��Ӧ����Ӧ����Ϣ
							for (int is = 0; is < IP_COUNT_1D; ++is)
							{
								element.IP(is).StrainLaststep() = element.IP(is).StrainCurrent();
								element.IP(is).StressLaststep() = element.IP(is).StressCurrent();
							}
						});

					//	����FEM�ڵ���Ϣ
					parallel_for_each(nids.begin(), nids.end(), [&](int nid)
						{
							TPdNode& node = pdModel.PdMeshCore().Node(nid);
							/************************************************************************/
							/* ��ǰ��������λ���ۼӵ�ȫ�ֽڵ���λ��                                 */
							/* �ڵ㵱ǰ����λ�Ƹ���Ϊ��һ��λ��+��ǰ��������λ��                    */
							/************************************************************************/
							node.TotalDisplacement() += node.IncrementalDisplacement();
							node.Coordinate() += node.IncrementalDisplacement().block(0, 0, 3, 1);
						});
				}
				/************************************************************************/
				/* End of FEM															*/
				/************************************************************************/


				/************************************************************************/
				/* Begin of PD&FEM														*/
				/************************************************************************/
				void				genGlobalMass()
				{
					double start = clock();
					TPdModel& pdModel = *m_pPdModel;

					const set<int>& nids = pdModel.PdMeshCore().GetNodeIdsByAll();
					int nCount = (int)(nids.size());
					m_GM.clear();
					m_GM.resize(DOF * nCount);

					for (const TPart& part : pdModel.Parts())
					{
						const set<int> eids = part.GetElementIds();
						const TMaterial& material = pdModel.Material(part.MaterialId());
						const TSection& section = pdModel.Section(part.SectionId());
			
						double Rho = material.GetMatValue("Rho");	
						double Iyy = section.GetSectionValue("ISS");
						double Izz = section.GetSectionValue("ITT");
						double Ixx = Iyy + Izz;

						/************************************************************************/
						/* ��װ��������                                                         */
						/************************************************************************/
						for (int ei : eids)
						{
							const TPdElement& element = pdModel.PdMeshCore().Element(ei);
							//	����ԪI�͵�ԪJ��Ӧ�Ľڵ����һ��nids������Ϊ8
							const vector<int>& nids = element.NodeIds();
							double A = element.Area();
							double L = element.SideLength();
							Eigen::Matrix<double, 6, 6> LocalM;
							LocalM.setZero();
							LocalM(0, 0) = Rho * A * L / 2.0;
							LocalM(1, 1) = Rho * A * L / 2.0;
							LocalM(2, 2) = Rho * A * L / 2.0;
							LocalM(3, 3) = Rho * Ixx * L / 2.0;
							LocalM(4, 4) = Rho * Iyy * L / 2.0;
							LocalM(5, 5) = Rho * Izz * L / 2.0;

							Eigen::Matrix<double, 6, 6> T;
							T.setZero();
							T.block(0, 0, 3, 3) = element.LocalCoorSystem();
							T.block(3, 3, 3, 3) = element.LocalCoorSystem();

							Eigen::Matrix<double, 6, 6> GlobalM = T * LocalM * T.transpose();

							for (int ni : nids)
							{
								int CurPos = ni * DOF;

								for (int col = 0; col < 6; ++col)
								{
									m_GM[CurPos][CurPos] += GlobalM(0, col);
									m_GM[CurPos + 1][CurPos + 1] += GlobalM(1, col);
									m_GM[CurPos + 2][CurPos + 2] += GlobalM(2, col);

									m_GM[CurPos + 3][CurPos + 3] += GlobalM(3, col);
									m_GM[CurPos + 4][CurPos + 4] += GlobalM(4, col);
									m_GM[CurPos + 5][CurPos + 5] += GlobalM(5, col);
								}
							}
						}
					}
					double total_time = (clock() - start) / 1000;
					cout << "genGlobalMassPD():\t\t" << total_time << endl;
				}
				bool				isConvergenced()
				{
					bool CONVERGENCED = true;

					TPdModel& pdModel = *m_pPdModel;
					set<int> nids = pdModel.PdMeshCore().GetNodeIdsByAll();

					double RatioIterIncr = -1E20;
					
					for (int nid : nids)
					{
						TPdNode& node = pdModel.PdMeshCore().Node(nid);
						const vector<TBoundaryPrescribedMotion>& BPMs = node.BoundaryPreMotion();
						//	ʩ��ǿ��λ�ƵĽڵ㲻��Ϊƽ���жϵĵ�
						if (BPMs.size() == 0)
						{
							double IterDisp = Module(node.IteratorDisplacement());
							double IncrDisp = Module(node.IncrementalDisplacement());
							double ratio_iter_incr = IterDisp / IncrDisp;
							if (ratio_iter_incr > RatioIterIncr)
							{
								RatioIterIncr = ratio_iter_incr;
							}
						}
					}

					if (RatioIterIncr > CONVERGENCE_FACTOR)
					{
						CONVERGENCED = false;
					}
					else
					{
						CONVERGENCED = true;
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
					Res.resize(4, 12);
					Res.setZero(); 
					 
					Res(0, 0) = -1.0 / L;
					Res(0, 6) = 1.0 / L;

					Res(1, 2) = (12 * chi - 6) / (L * L);
					Res(1, 4) = -(6 * chi - 4) / L;
					Res(1, 8) = -Res(1, 2);
					Res(1, 10) = -(6 * chi - 2) / L;

					Res(2, 1) = (12 * chi - 6) / (L * L);
					Res(2, 5) = (6 * chi - 4) / L;
					Res(2, 7) = -Res(2, 1);
					Res(2, 11) = (6 * chi - 2) / L;

					Res(3, 3) = -1.0 / L;
					Res(3, 9) = 1.0 / L;

					return Res;
				}
				Eigen::MatrixXd		BN_star_Matrix(double L, double chi, const MatrixXd& delta_u)
				{
					Eigen::MatrixXd Res;
					Res.resize(4, 12);
					Res.setZero();
					Eigen::MatrixXd G = G_Matrix(L, chi);

					Eigen::MatrixXd delta_A;
					delta_A.resize(4, 2);
					delta_A.setZero();
					delta_A.block(0, 0, 1, 2)  = (G * delta_u).transpose(); 

					Res = 0.5 * delta_A * G;
					
					return Res;
				}
				/************************************************************************/
				/* End of PD&FEM														*/
				/************************************************************************/							
			private:
				vector< map<int, double> > m_GK;			//	Stiffness matrix
				vector< map<int, double> > m_GM;			//	Mass matrix	

				Eigen::Matrix4d				D_Elastic;		//	Elastic matrix
			private:
				TPdModel*					m_pPdModel;
			};
		} // END OF PERIDYNAMIC
	} // END OF DLUT
} // END OF SAE

#endif