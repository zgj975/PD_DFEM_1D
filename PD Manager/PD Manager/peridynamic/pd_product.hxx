/*==============================================================================

Copyright 2017 Dalian University of Technology .
All rights reserved

================================================================================
-- Please append file description informations here --
PdProduct
1. Set calculate parameters for SOLVER
2. TPdNode can be : TPdNodeBBPD / TPdNodeNOSPD / TPdNodeBeamPD
3. TSolver can be : TSolverBBPD / TSolverNOSPD / TSolverBeamPD
================================================================================
Date            Name                    Description of Change
2017/12/24		Zheng Guojun			Create
2019/07/31		Zheng Guojun			Unified the solver call process
$HISTORY$
================================================================================
*/

#ifndef DLUT_SAE_PERIDYNAMIC_PD_PRODUCT_HXX_20171224
#define DLUT_SAE_PERIDYNAMIC_PD_PRODUCT_HXX_20171224

#include "pd_serialize.hxx"
#include "pd_solver_beampd.hxx"

#include <time.h>

namespace DLUT
{
	namespace SAE
	{
		namespace PERIDYNAMIC
		{
			template<typename TSolver>
			class PdProduct
			{
			public:
				PdProduct() {}
				~PdProduct() {}

			public:
				void		SetFilePath(string ip_lsdyna_info)
				{
					m_str_lsdyna_file = ip_lsdyna_info;
					string::size_type index = m_str_lsdyna_file.find_last_of('.');
					if (m_str_lsdyna_file.npos != index)
					{
						m_str_output_ascii = m_str_lsdyna_file.substr(0, index) + "_m_" + to_string(RATIO_OF_HORIZON_MESHSIZE) + ".ascii";
						m_str_output_nodes = m_str_lsdyna_file.substr(0, index) + "_nodout.dat";
						m_str_output_elems = m_str_lsdyna_file.substr(0, index) + "_eleout.dat";
						m_str_modal_ascii = m_str_lsdyna_file.substr(0, index) + ".modal";
					}
				}
				void		SetRatioOfHorizonMeshsize(double ratio_h_m)
				{
					DLUT::SAE::PERIDYNAMIC::RATIO_OF_HORIZON_MESHSIZE = ratio_h_m;
				}
				void		SetConstantHorizon(double horizon)
				{
					DLUT::SAE::PERIDYNAMIC::HORIZON = horizon;
				}
				void		SetUsedConstantHorizon(bool used)
				{
					DLUT::SAE::PERIDYNAMIC::USE_CONSTANT_HORIZON = used;
				}
				void		SetLoadStep(int load_step)
				{
					m_load_step = load_step;
					if (m_load_step < 1)
					{
						m_load_step = 1;
					}
				}
				void		SetMaxIteratorNums(int max_iter_nums)
				{
					DLUT::SAE::PERIDYNAMIC::MAX_INTERATOR_NUMS = max_iter_nums;
				}
				void		SetConvergenceFactor(double conv_factor)
				{
					DLUT::SAE::PERIDYNAMIC::CONVERGENCE_FACTOR = conv_factor;
				}
				void		SetTimeStep(double time_step)
				{
					m_time_step = time_step;
					if (m_time_step < 0)
					{
						m_time_step = 1;
					}
				}
				void		SetPlotFrames(int plot_frames)
				{
					m_plot_frames = plot_frames;
				}
				void		SetIteratorNums(int iterator_nums)
				{								
					m_iterator_nums = iterator_nums;					
				}
			public:
				void		ExplicitAnalysis()
				{		
					cout << "************************************************************************" << endl;
					cout << "\t\tEXPLICIT ANALYSIS " << endl;
					cout << "************************************************************************" << endl;
					double start = clock();

					int plot_frames = m_iterator_nums / m_plot_frames;
					//	Step 1: Read Pre-Process information
					fun_initPreData();
					seri_ascii.OutputHwAscii(0);
						
					//	Step 2: Refresh initial velocity/Acceleration/
					m_pd_model.RefreshInitVelocityInfo();

					TSolver solver;
					solver.Attach(m_pd_model);
					//	Step 4: Solution						
					cout << "************************************************************************" << endl;
					double cur_step_start = clock();
					for (int tt = 1; tt <= m_iterator_nums; ++tt)
					{
						//	Step 3: Refresh body force
						m_pd_model.RefreshLoadNodeInfo(tt * m_time_step);

						//	Step 4.0: Load
						m_pd_model.RefreshPreMotionInfo(tt, m_time_step);

						//	Step 4.1: Calculate all the PD force 
						solver.ExplicitSolve(m_time_step);

						//	Step 4.2 刷新约束信息
						m_pd_model.RefreshBoundarySpcInfo();

						//	Step 4.4: Check fracture
						solver.RefreshFracture();

						//	Step 4.5: Output information (INTERVAL = tt*m_dt)
						if ((tt % plot_frames) == 0)
						{
							seri_ascii.OutputHwAscii(tt);
							
							double cur_step_cost = (clock() - cur_step_start) / 1000;
							cur_step_start = clock();
							cout << "ITERATOR NUMS = " << setw(5) << tt << ", ";
							cout << "CALCULATE TIME = " << setw(5) << cur_step_cost << endl;
						}				

						seri_node.OutputNodesInfo(tt);
						seri_element.OutputElementInfo(tt);
					}
					
					fun_finishOut();
					double total_time = (clock() - start) / 1000;
					cout << "************************************************************************" << endl;
					cout << "COMPLETED, TOTAL TIME = " << total_time << endl;
				}
				void		ImplicitAnalysis()
				{
					cout << "************************************************************************" << endl;
					cout << "\t\tIMPLICIT ANALYSIS " << endl;
					cout << "************************************************************************" << endl;
					double start = clock();

					//	Step 1: Read Pre-Process information
					fun_initPreData();

					seri_ascii.OutputHwAscii(0);

					TSolver solver;
					solver.Attach(m_pd_model);
					
					cout << "************************************************************************" << endl;
					for (int cur_step = 0; cur_step <= m_load_step; ++cur_step)
					{
						double cur_step_start = clock();

						int ITERATOR_NUM = 0;
						bool is_convergenced = solver.ImplicitSolve(cur_step, m_load_step, ITERATOR_NUM);
						double cur_step_cost = (clock() - cur_step_start) / 1000;

						seri_ascii.OutputHwAscii(cur_step);

						cout << "CURRENT STEP = " << setw(5) << cur_step << ", ";
						cout << "ITERATOR NUMS = " << setw(5) << ITERATOR_NUM << ", ";
						cout << "CALCULATE TIME = " << setw(5) << cur_step_cost << ", ";
						if (is_convergenced)
						{
							cout << "Convergenced.\n";
						}
						else
						{
							cout << "Not convergenced.\n";

							break;
						}
					}

					fun_finishOut();

					double total_time = (clock() - start) / 1000;
					cout << "************************************************************************" << endl;
					cout << "COMPLETED, TOTAL TIME = " << total_time << endl;
				}
			private:
				void		fun_initPreData()
				{
//					cout << "Begin to read Pre-Process datas..." << endl;
					//	Read information from *KEY files (LSDYNA)
					m_pd_model.Initialize();
					SerializeLsdyna dyna_seri;
					dyna_seri.ReadLsdynaFile(m_str_lsdyna_file, m_pd_model);

					//	读入信息之后更新PART的信息
 					m_pd_model.UpdatePartInfo();
					//	Update family information for all PD nodes
					m_pd_model.UpdateFamilyInParts();

					//	Update initial damage value for all PD nodes
					m_pd_model.GenerateInitDamage();
					//	Generate initial crevice and delete the corresponding family information
					m_pd_model.GenerateInitCrevice();
//					cout << "Finished to read Pre-Process datas..." << endl;

					//	Initialize output serialize.
					seri_ascii.AttachModel(m_pd_model);
					seri_ascii.AttachFile(m_str_output_ascii);

					seri_element.AttachModel(m_pd_model);
					seri_element.AttachFile(m_str_output_elems);

					seri_node.AttachModel(m_pd_model);
					seri_node.AttachFile(m_str_output_nodes);

					string::size_type index = m_str_lsdyna_file.find_last_of('.');
					if (m_str_lsdyna_file.npos != index)
					{
						string lsdyna_file_new = m_str_lsdyna_file.substr(0, index) + "_discret.key";
						dyna_seri.WriteLsdynaFile(lsdyna_file_new, m_pd_model);
					}
				}
				void		fun_finishOut()
				{
					seri_ascii.DetachFile();
					seri_element.DetachFile();
					seri_node.DetachFile();					
				}
			private:
				TPdModel	m_pd_model;				//	PD data of everything
			private:
				int			m_load_step;			//	Numbers of Load step
				double		m_time_step;			//	Time step
				int			m_iterator_nums;		//	Numbers of Iterator at every load step
				int			m_plot_frames;			//	Numbers of time interval for output

			private:
				string		m_str_lsdyna_file;		//	Filename of other information
				string		m_str_output_ascii;		//	Filename of output information
				string		m_str_modal_ascii;		//	Filename of modal information
				string		m_str_output_nodes;		//	Filename of output NODES information
				string		m_str_output_elems;		//	Filename of output ELEMS information
			private:
				SerializeHwAscii		seri_ascii;
				SerializeElementsInfo	seri_element;
				SerializeNodesInfo		seri_node;
			};
		} // END OF PERIDYNAMIC
	} // END OF DLUT
} // END OF SAE

#endif