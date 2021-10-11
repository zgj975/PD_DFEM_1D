#ifndef DLUT_SAE_PERIDYNAMIC_FEM_SHAPE_FUNCTIONS_20210528
#define DLUT_SAE_PERIDYNAMIC_FEM_SHAPE_FUNCTIONS_20210528

namespace DLUT
{
	namespace SAE
	{
		namespace PERIDYNAMIC
		{
			//	Beam单元形函数，其中Chi范围为[0,1]
			Eigen::MatrixXd		N_SF_BEAM(double L, double Chi)
			{
				Eigen::MatrixXd Res;
				Res.resize(6, 12);
				Res.setZero();
				Res(0, 0) = 1 - Chi;
				Res(0, 6) = Chi;

				Res(1, 1) = 1 - 3 * Chi * Chi + 2 * Chi * Chi * Chi;
				Res(1, 5) = (Chi - 2 * Chi * Chi + Chi * Chi * Chi) * L;
				Res(1, 7) = 3 * Chi * Chi - 2 * Chi * Chi * Chi;
				Res(1, 11) = (Chi * Chi * Chi - Chi * Chi) * L;

				Res(2, 2) = 1 - 3 * Chi * Chi + 2 * Chi * Chi * Chi;
				Res(2, 4) = -(Chi - 2 * Chi * Chi + Chi * Chi * Chi) * L;
				Res(2, 8) = 3 * Chi * Chi - 2 * Chi * Chi * Chi;
				Res(2, 10) = -(Chi * Chi * Chi - Chi * Chi) * L;

				Res(3, 3) = 1 - Chi;
				Res(3, 9) = Chi;

				Res(4, 2) = -(-6 * Chi + 6 * Chi * Chi) / L;
				Res(4, 4) = (1 - 4 * Chi + 3 * Chi * Chi);
				Res(4, 8) = -(6 * Chi - 6 * Chi * Chi) / L;
				Res(4, 10) = (3 * Chi * Chi - 2 * Chi);

				Res(5, 1) = (-6 * Chi + 6 * Chi * Chi) / L;
				Res(5, 5) = (1 - 4 * Chi + 3 * Chi * Chi);
				Res(5, 7) = (6 * Chi - 6 * Chi * Chi) / L;
				Res(5, 11) = (3 * Chi * Chi - 2 * Chi);

				return Res;
			}

			//	Rod单元形函数，其中Chi范围为[0,1]
			Eigen::MatrixXd		N_SF_ROD(double Chi)
			{
				Eigen::MatrixXd Res;
				Res.resize(3, 6);
				Res.setZero();
				Res(0, 0) = 1 - Chi;
				Res(0, 3) = Chi;
				Res(1, 1) = 1 - Chi;
				Res(1, 4) = Chi;
				Res(2, 2) = 1 - Chi;
				Res(2, 5) = Chi;

				return Res;
			}
		}
	}
}

#endif