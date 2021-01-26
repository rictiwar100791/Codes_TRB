//Reference:  Tiwari, Jayaswal, Sinha (2021). Competitive Hub Location Problem: Model and Solution Approaches, Transportation Research-B, Under Review.
//This code uses three methods: Mixed Integer Second Order Conic Program (MISOCP), Cutting Plane Algorithm (CPA) and the Lagrangian Relaxation with Cutting Plane Algorithm (LR-CPA).
//This code corresponds to the multiple allocation (MA) version of the problem as described in Section 3.2. Also, the experimenst are run on standard data-sets of CAB 10, 15, 20 and 25, AP50.


#pragma warning(disable : 4996)//For Visual Studio 2012
#include<stdio.h>
#include<conio.h>
#include<iostream>
#include<fstream>
#include<iosfwd>
#include<string>
#include <deque>
#include <sstream>
#include <time.h>
#include <stdlib.h>
#define TIMLIM 3600.0

#include <ilcplex/ilocplex.h>
#include <ilconcert/ilosys.h>

ILOSTLBEGIN
IloEnv env;
typedef IloArray<IloNumArray> Num2DMatrix;
typedef IloArray<Num2DMatrix>Num3DMatrix;
typedef IloArray<Num3DMatrix>Num4DMatrix;
typedef IloArray<IloBoolVarArray> BoolVar2DMatrix;
typedef IloArray<IloBoolArray> Bool2DMatrix;
typedef IloArray<IloNumVarArray> NumVar2DMatrix;
typedef IloArray<NumVar2DMatrix> NumVar3DMatrix;
typedef IloArray<NumVar3DMatrix> NumVar4DMatrix;

typedef IloArray<IloNumArray> Num2DMatrix;
IloInt N,c;
Num2DMatrix h(env);//flow from i to j nodes.
Num2DMatrix dist_xy(env);
typedef struct XVAL 
{
	int index;
	double value;
} XVAL;

int Comparevalue(const void *a, const void *b)
{
	if (((XVAL *)a)->value<((XVAL *)b)->value)
		return 1;
	if (((XVAL *)a)->value>((XVAL *)b)->value)
		return -1;
	return 0;
}

void SOCP(int p, double al) //Refer to Section 4, Mixed Integer Second Order Conic Program (MISOCP)
{
	IloInt A;//Number of nodes and hubs
	IloNum eps;//later assigned as eps = cplex.getParam(IloCplex::EpInt); EpInt is the tolerance gap for integer variables
	IloNum be, ga, de;
	A = 1000;
	be = de = 1;
	ga = .75;
	IloNum c_collect, c_distribute;
	c_collect = c_distribute = 1;
	IloModel model(env);
	IloCplex cplex(model);
	Num4DMatrix cost_xy(env, N);
	Num4DMatrix time_xy(env, N);
	Num4DMatrix s(env, N);
	Num4DMatrix o(env, N);
	clock_t start;

	////////////////////////DEFINE VARIABLES///////////////////////

	for (int i = 0; i < N; i++)
	{
		cost_xy[i] = Num3DMatrix(env, N);
		time_xy[i] = Num3DMatrix(env, N);
		s[i] = Num3DMatrix(env, N);
		o[i] = Num3DMatrix(env, N);
		for (int j = 0; j < N; j++)
		{
			s[i][j] = Num2DMatrix(env, N);
			o[i][j] = Num2DMatrix(env, N);
			cost_xy[i][j] = Num2DMatrix(env, N);
			time_xy[i][j] = Num2DMatrix(env, N);

			for (int k = 0; k < N; k++)
			{
				s[i][j][k] = IloNumArray(env, N);
				o[i][j][k] = IloNumArray(env, N);
				cost_xy[i][j][k] = IloNumArray(env, N);
				time_xy[i][j][k] = IloNumArray(env, N);
			}
		}
	}


	cout << "Attractiveness index is " << A << endl;
	cout << "cost decay parameter is " << be << endl;
	cout << "time decay parameter is " << de << endl;
	eps = cplex.getParam(IloCplex::EpInt);

	for (int i = 0; i<N; i++)
	{
		for (int j = 0; j<N; j++)
		{
			for (int k = 0; k<N; k++)
			{
				for (int l = 0; l<N; l++)
				{
					if (i != j)
					{
						cost_xy[i][j][k][l] = c_collect*dist_xy[i][k] + al*dist_xy[k][l] + c_distribute*dist_xy[l][j];
						time_xy[i][j][k][l] = dist_xy[i][k] + dist_xy[k][l] + dist_xy[l][j];
						if (k != l)
						{
							s[i][j][k][l] = A / (ga*pow(time_xy[i][j][k][l], be) + (1 - ga)*pow(cost_xy[i][j][k][l], de));
							o[i][j][k][l] = A / (ga*pow(time_xy[i][j][k][l], be) + (1 - ga)*pow(cost_xy[i][j][k][l], de));
						}
						else if (k == l)
						{
							s[i][j][k][l] = A*1.25 / (ga*pow(time_xy[i][j][k][l], be) + (1 - ga)*pow(cost_xy[i][j][k][l], de));
							o[i][j][k][l] = A*1.25 / (ga*pow(time_xy[i][j][k][l], be) + (1 - ga)*pow(cost_xy[i][j][k][l], de));
						}
					}
					if (i == j)
					{
						cost_xy[i][j][k][l] = 0;
						time_xy[i][j][k][l] = 0;
						s[i][j][k][l] = eps;
						o[i][j][k][l] = eps;
					}
				}
			}
		}
	}

	//Define decision variable
	NumVar4DMatrix X(env, N);//1, if hubs are located at both k and l
	Num4DMatrix X_val(env, N);//1, if hubs are located at both k and l
	IloNumArray Y_val(env, N);
	IloBoolVarArray Y(env, N);
	IloNumArray best_Y(env, N);
	IloBoolVarArray comp_Y(env, N);
	IloNumArray compY_val(env, N);
	NumVar4DMatrix comp_X(env, N);//1, if hubs are located at both k and l
	Num4DMatrix comp_X_val(env, N);//1, if hubs are located at both k and l
	NumVar2DMatrix Z(env, N);
	Num2DMatrix Z_val(env, N);
	NumVar2DMatrix R(env, N);
	Num2DMatrix R_val(env, N);
	IloNumVarArray t(env, 3, 0, IloInfinity, ILOFLOAT);
	IloNumArray t_val(env, 3);
	IloNumVar n(env, 0, IloInfinity, ILOFLOAT);
	Num3DMatrix f(env, N);
	NumVar2DMatrix r(env, N);
	Num2DMatrix r_val(env, N);

	//Assigning value to variables

	for (int i = 0; i < N; i++)
	{
		X[i] = NumVar3DMatrix(env, N);
		X_val[i] = Num3DMatrix(env, N);
		Z[i] = IloNumVarArray(env, N, 0, 1, ILOFLOAT);
		Z_val[i] = IloNumArray(env, N);
		R[i] = IloNumVarArray(env, N, 0, IloInfinity, ILOFLOAT);
		R_val[i] = IloNumArray(env, N);
		f[i] = Num2DMatrix(env, N);
		r[i] = IloNumVarArray(env, N, 0, IloInfinity, ILOFLOAT);
		r_val[i] = IloNumArray(env, N);
		for (int j = 0; j < N; j++)
		{
			X[i][j] = NumVar2DMatrix(env, N);
			X_val[i][j] = Num2DMatrix(env, N);
			f[i][j] = IloNumArray(env, 3);
			for (int k = 0; k < N; k++)
			{
				X[i][j][k] = IloNumVarArray(env, N, 0, IloInfinity, ILOFLOAT);
				X_val[i][j][k] = IloNumArray(env, N);
			}
		}
	}



	for (int i = 0; i < N; i++)
	{
		comp_X[i] = NumVar3DMatrix(env, N);
		comp_X_val[i] = Num3DMatrix(env, N);
		for (int j = 0; j < N; j++)
		{
			comp_X[i][j] = NumVar2DMatrix(env, N);
			comp_X_val[i][j] = Num2DMatrix(env, N);
			for (int k = 0; k < N; k++)
			{
				comp_X[i][j][k] = IloNumVarArray(env, N, 0, 1, ILOFLOAT);
				comp_X_val[i][j][k] = IloNumArray(env, N);
			}
		}
	}

	IloModel model_comp(env);
	IloCplex cplex_comp(model_comp);
	IloExpr Obj_comp(env); // Creates an expression with the name Obj (Objective)

	for (int i = 0; i<N; i++)

	{

		for (int j = 0; j<N; j++)

		{

			for (int k = 0; k<N; k++)

			{

				for (int m = 0; m<N; m++)

				{

					Obj_comp += (c_collect*dist_xy[i][k] + al*dist_xy[k][m] + c_distribute*dist_xy[m][j]) *comp_X[i][j][k][m];

				}

			}

		}

	}


	model_comp.add(IloMinimize(env, Obj_comp)); // IloMinimize is used for minimization problems

	Obj_comp.end(); // Clear out the memory allocated for the expression 



	//Constraint 1: 

	for (int i = 0; i<N; i++)

	{

		for (int j = 0; j<N; j++)

		{
			IloExpr SumX(env);

			for (int m = 0; m < N; m++)

			{

				for (int k = 0; k < N; k++)

				{
					SumX += comp_X[i][j][k][m];
				}

			}
			model_comp.add(SumX == 1);
			SumX.end();

		}

	}

	//Constraint 2: sum {k in 1..N}(Y[k])=p;

	IloExpr SumZkk(env);

	for (int k = 0; k<N; k++)

	{

		SumZkk += comp_Y[k];

	}

	model_comp.add(SumZkk == p);

	SumZkk.end();



	//Constraint 3: for Facet Inducing

	for (int i = 0; i<N; i++)

	{

		for (int j = 0; j<N; j++)

		{

			for (int k = 0; k<N; k++)

			{

				IloExpr SumX(env);

				for (int m = 0; m < N; m++)

				{
					SumX += comp_X[i][j][k][m];

					if (m != k)
					{
						SumX += comp_X[i][j][m][k];
					}
				}

				model_comp.add(SumX <= comp_Y[k]);
				SumX.end();

			}

		}

	}

	//==============================================================================

	//Optimize

	cplex_comp.setOut(env.getNullStream()); //This is to supress the output of Branch & Bound Tree on screen

	cplex_comp.setWarning(env.getNullStream()); //This is to supress warning messages on screen

	cplex_comp.solve();//solving the MODEL

	if (cplex_comp.getStatus() == IloAlgorithm::Infeasible) // if the problem is infeasible

	{

		cout << "Problem is Infeasible" << endl;

	}

	for (int k = 0; k<N; k++)

	{

		if (cplex_comp.getValue(comp_Y[k]) > 1 - eps)

		{

			compY_val[k] = 1;

		}

		else

		{

			compY_val[k] = 0;

		}

	}

	for (int i = 0; i<N; i++)

	{

		for (int j = 0; j<N; j++)

		{

			for (int k = 0; k<N; k++)

			{

				for (int m = 0; m<N; m++)

				{

					if (cplex_comp.getValue(comp_X[i][j][k][m]) > 1 - eps)

					{

						comp_X_val[i][j][k][m] = 1;

					}

					else

					{

						comp_X_val[i][j][k][m] = 0;

					}

				}

			}

		}

	}

	cout << "Open Hubs for the competitor are==========================" << endl;
	for (int k = 0; k < N; k++)
	{
		if (compY_val[k] > 0)
		{
			cout << k + 1 << endl;
		}
	}

	//====================================Entrant's problem================================================================//


	/*
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			for (int g = 0; g < N; g++)
			{
				if (g == 0)
				{
					f[i][j][g] = h[i][j] * .3;
				}
				if (g == 1)
				{
					f[i][j][g] = h[i][j] * .6;
				}
				if (g == 2)
				{
					f[i][j][g] = h[i][j] * .9;
				}

			}
		}
	}
	*/

	IloNum UB = IloInfinity;

	//Begin Iterations==============================================================================

	start = clock();

	//IloNumVar theta(env, 0, IloInfinity, ILOFLOAT);
	IloExpr Obj(env); // Creates an expression with the name Obj (Objective)
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			Obj += h[i][j] * Z[i][j];
		}
	}
	
	//Obj -= omega*n;

	// model.add is a function to add constraints and objectives in the CPLEX environment
	model.add(IloMaximize(env, Obj)); // IloMinimize is used for minimization problems
	Obj.end(); // Clear out the memory allocated for the expression 
	// or is ||

	IloExpr SumYp(env);

	for (int k = 0; k<N; k++)

	{

		SumYp += Y[k];

	}

	model.add(SumYp == p);

	SumYp.end();

	//Constraint 2:
	for (int i = 0; i<N; i++)

	{

		for (int j = 0; j<N; j++)

		{
			IloExpr SumX(env);
			for (int m = 0; m < N; m++)

			{

				for (int k = 0; k < N; k++)

				{

					SumX += X[i][j][k][m];


				}

			}
			model.add(SumX == 1);
			SumX.end();

		}

	}




	//Constraint 3: for Facet Inducing

	for (int i = 0; i<N; i++)

	{

		for (int j = 0; j<N; j++)

		{

			for (int k = 0; k<N; k++)

			{

				IloExpr SumX(env);

				for (int m = 0; m < N; m++)

				{
					SumX += X[i][j][k][m];

					if (m != k)
					{
						SumX += X[i][j][m][k];
					}
				}

				model.add(SumX <= Y[k]);
				SumX.end();

			}

		}

	}



	//Constraint 5


	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			IloExpr Sum_o(env);
			for (int m = 0; m<N; m++)
			{
				for (int n = 0; n<N; n++)
				{
					Sum_o += X[i][j][m][n] * o[i][j][m][n] + s[i][j][m][n] * comp_X_val[i][j][m][n];
				}
			}
			model.add(r[i][j] == Sum_o);
			Sum_o.end();
		}
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			IloExpr Sum_o(env);
			for (int m = 0; m<N; m++)
			{
				for (int n = 0; n<N; n++)
				{
					Sum_o += s[i][j][m][n] * comp_X_val[i][j][m][n] * R[i][j];
				}
			}
			model.add(Z[i][j] + Sum_o == 1);
			Sum_o.end();
		}
	}

	/*
	for (int g = 0; g < 3; g++)
	{
		IloExpr Sumt(env);
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				Sumt += f[i][j][g] * Z[i][j];
			}
		}

		model.add(Sumt >= t[g]);
		Sumt.end();
	}

	IloExpr Sumn(env);
	for (int g = 0; g < 3; g++)
	{
		Sumn += t[g] * t[g];
	}
	model.add(Sumn <= n*n);
	Sumn.end();

	*/

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			IloExpr Sum_o(env);
			for (int m = 0; m<N; m++)
			{
				for (int n = 0; n<N; n++)
				{
					Sum_o += o[i][j][m][n] * X[i][j][m][n] * X[i][j][m][n];
				}
			}
			model.add((Z[i][j] + r[i][j])*(Z[i][j] + r[i][j]) >= 2 * Sum_o + Z[i][j] * Z[i][j] + r[i][j] * r[i][j]);
			Sum_o.end();
		}
	}


	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{

			model.add((R[i][j] + r[i][j])*(R[i][j] + r[i][j]) >= 2 + R[i][j] * R[i][j] + r[i][j] * r[i][j]);

		}
	}

	cplex.setParam(IloCplex::TiLim, TIMLIM);
	// Print results
	IloCplex cplex_(model);
	//cplex.setOut(env.getNullStream()); // if we get an empty stream in return
	//cplex.setWarning(env.getNullStream());
	cplex.solve();//solving the MODEL


	if (cplex.getStatus() == IloAlgorithm::Infeasible) // if the problem is infeasible
	{
		cout << "Problem Infeasible" << endl;
	}
	UB = cplex.getObjValue();

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			for (int m = 0; m < N; m++)
			{
				for (int n = 0; n < N; n++)
				{
					X_val[i][j][m][n] = cplex.getValue(X[i][j][m][n]);
				}
			}
		}
	}

	for (int k = 0; k < N; k++)
	{
		Y_val[k] = cplex.getValue(Y[k]);
	}


	IloNum MS = 0;
	IloNum Flow = 0;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			Flow = Flow + h[i][j];
		}

	}
	MS = UB / Flow;
	cout << "==========================Original problem==================================" << endl;

	cout << "The objective function to the MACOHLP is" << UB << endl;
	cout << "The market share captured by the entrant is " << MS << endl;
	cout << "Total Run time \t" << (clock() - start) / (double)CLOCKS_PER_SEC << endl;

}




void CPA(int p, double al) //Refer to Appendix, Cutting Plane Algorithm (CPA)
{

	IloInt A;//Number of nodes and hubs
	IloNum eps;//later assigned as eps = cplex.getParam(IloCplex::EpInt); EpInt is the tolerance gap for integer variables
	IloNum be, ga, de;
	A = 1000;
	be = de = 1;
	ga = .75;
	IloNum c_collect, c_distribute;
	c_collect = c_distribute = 1;
	IloModel model(env);
	IloCplex cplex(model);
	Num4DMatrix cost_xy(env, N);
	Num4DMatrix time_xy(env, N);
	Num4DMatrix s(env, N);
	Num4DMatrix o(env, N);
	clock_t start, run_SB = 0, run_SB_tot = 0, start_iter;

	////////////////////////DEFINE VARIABLES///////////////////////

	for (int i = 0; i < N; i++)
	{
		cost_xy[i] = Num3DMatrix(env, N);
		time_xy[i] = Num3DMatrix(env, N);
		s[i] = Num3DMatrix(env, N);
		o[i] = Num3DMatrix(env, N);
		for (int j = 0; j < N; j++)
		{
			s[i][j] = Num2DMatrix(env, N);
			o[i][j] = Num2DMatrix(env, N);
			cost_xy[i][j] = Num2DMatrix(env, N);
			time_xy[i][j] = Num2DMatrix(env, N);

			for (int k = 0; k < N; k++)
			{
				s[i][j][k] = IloNumArray(env, N);
				o[i][j][k] = IloNumArray(env, N);
				cost_xy[i][j][k] = IloNumArray(env, N);
				time_xy[i][j][k] = IloNumArray(env, N);
			}
		}
	}


	cout << "Attractiveness index is " << A << endl;
	cout << "cost decay parameter is " << be << endl;
	cout << "time decay parameter is " << de << endl;
	eps = cplex.getParam(IloCplex::EpInt);

	for (int i = 0; i<N; i++)
	{
		for (int j = 0; j<N; j++)
		{
			for (int k = 0; k<N; k++)
			{
				for (int l = 0; l<N; l++)
				{
					if (i != j)
					{
						cost_xy[i][j][k][l] = c_collect*dist_xy[i][k] + al*dist_xy[k][l] + c_distribute*dist_xy[l][j];
						time_xy[i][j][k][l] = dist_xy[i][k] + dist_xy[k][l] + dist_xy[l][j];
						if (k != l)
						{
							s[i][j][k][l] = A / (ga*pow(time_xy[i][j][k][l], be) + (1 - ga)*pow(cost_xy[i][j][k][l], de));
							o[i][j][k][l] = A / (ga*pow(time_xy[i][j][k][l], be) + (1 - ga)*pow(cost_xy[i][j][k][l], de));
						}
						else if (k == l)
						{
							s[i][j][k][l] = A*1.25 / (ga*pow(time_xy[i][j][k][l], be) + (1 - ga)*pow(cost_xy[i][j][k][l], de));
							o[i][j][k][l] = A*1.25 / (ga*pow(time_xy[i][j][k][l], be) + (1 - ga)*pow(cost_xy[i][j][k][l], de));
						}
					}
					if (i == j)
					{
						cost_xy[i][j][k][l] = 0;
						time_xy[i][j][k][l] = 0;
						s[i][j][k][l] = eps;
						o[i][j][k][l] = eps;
					}
				}
			}
		}
	}

	//Define decision variable
	NumVar4DMatrix X(env, N);//1, if hubs are located at both k and l
	Num4DMatrix X_val(env, N);//1, if hubs are located at both k and l
	IloNumArray Y_val(env, N);
	IloBoolVarArray Y(env, N);
	IloNumArray best_Y(env, N);
	IloBoolVarArray comp_Y(env, N);
	IloNumArray compY_val(env, N);
	NumVar4DMatrix comp_X(env, N);//1, if hubs are located at both k and l
	Num4DMatrix comp_X_val(env, N);//1, if hubs are located at both k and l
	NumVar2DMatrix Z(env, N);
	Num2DMatrix Z_val(env, N);
	NumVar2DMatrix R(env, N);
	Num2DMatrix R_val(env, N);

	//Assigning value to variables

	for (int i = 0; i < N; i++)
	{
		X[i] = NumVar3DMatrix(env, N);
		X_val[i] = Num3DMatrix(env, N);
		Z[i] = IloNumVarArray(env, N, 0, 1, ILOFLOAT);
		Z_val[i] = IloNumArray(env, N);
		R[i] = IloNumVarArray(env, N, 0, IloInfinity, ILOFLOAT);
		R_val[i] = IloNumArray(env, N);
		for (int j = 0; j < N; j++)
		{
			X[i][j] = NumVar2DMatrix(env, N);
			X_val[i][j] = Num2DMatrix(env, N);
			for (int k = 0; k < N; k++)
			{
				X[i][j][k] = IloNumVarArray(env, N, 0, IloInfinity, ILOFLOAT);
				X_val[i][j][k] = IloNumArray(env, N);
			}
		}
	}

	for (int i = 0; i < N; i++)
	{
		comp_X[i] = NumVar3DMatrix(env, N);
		comp_X_val[i] = Num3DMatrix(env, N);
		for (int j = 0; j < N; j++)
		{
			comp_X[i][j] = NumVar2DMatrix(env, N);
			comp_X_val[i][j] = Num2DMatrix(env, N);
			for (int k = 0; k < N; k++)
			{
				comp_X[i][j][k] = IloNumVarArray(env, N, 0, 1, ILOFLOAT);
				comp_X_val[i][j][k] = IloNumArray(env, N);
			}
		}
	}

	IloModel model_comp(env);
	IloCplex cplex_comp(model_comp);
	IloExpr Obj_comp(env); // Creates an expression with the name Obj (Objective)

	for (int i = 0; i<N; i++)

	{

		for (int j = 0; j<N; j++)

		{

			for (int k = 0; k<N; k++)

			{

				for (int m = 0; m<N; m++)

				{

					Obj_comp += (c_collect*dist_xy[i][k] + al*dist_xy[k][m] + c_distribute*dist_xy[m][j]) *comp_X[i][j][k][m];

				}

			}

		}

	}


	model_comp.add(IloMinimize(env, Obj_comp)); // IloMinimize is used for minimization problems

	Obj_comp.end(); // Clear out the memory allocated for the expression 



	//Constraint 1: 

	for (int i = 0; i<N; i++)

	{

		for (int j = 0; j<N; j++)

		{
			IloExpr SumX(env);

			for (int m = 0; m < N; m++)

			{

				for (int k = 0; k < N; k++)

				{
					SumX += comp_X[i][j][k][m];
				}

			}
			model_comp.add(SumX == 1);
			SumX.end();

		}

	}

	//Constraint 2: sum {k in 1..N}(Y[k])=p;

	IloExpr SumZkk(env);

	for (int k = 0; k<N; k++)

	{

		SumZkk += comp_Y[k];

	}

	model_comp.add(SumZkk == p);

	SumZkk.end();



	//Constraint 3: for Facet Inducing

	for (int i = 0; i<N; i++)

	{

		for (int j = 0; j<N; j++)

		{

			for (int k = 0; k<N; k++)

			{

				IloExpr SumX(env);

				for (int m = 0; m < N; m++)

				{
					SumX += comp_X[i][j][k][m];

					if (m != k)
					{
						SumX += comp_X[i][j][m][k];
					}
				}

				model_comp.add(SumX <= comp_Y[k]);
				SumX.end();

			}

		}

	}

	//==============================================================================



	//Optimize

	cplex_comp.setOut(env.getNullStream()); //This is to supress the output of Branch & Bound Tree on screen

	cplex_comp.setWarning(env.getNullStream()); //This is to supress warning messages on screen

	cplex_comp.solve();//solving the MODEL

	if (cplex_comp.getStatus() == IloAlgorithm::Infeasible) // if the problem is infeasible

	{

		cout << "Problem is Infeasible" << endl;

	}

	for (int k = 0; k<N; k++)

	{

		if (cplex_comp.getValue(comp_Y[k]) > 1 - eps)

		{

			compY_val[k] = 1;

		}

		else

		{

			compY_val[k] = 0;

		}

	}

	for (int i = 0; i<N; i++)

	{

		for (int j = 0; j<N; j++)

		{

			for (int k = 0; k<N; k++)

			{

				for (int m = 0; m<N; m++)

				{

					if (cplex_comp.getValue(comp_X[i][j][k][m]) > 1 - eps)

					{

						comp_X_val[i][j][k][m] = 1;

					}

					else

					{

						comp_X_val[i][j][k][m] = 0;

					}

				}

			}

		}

	}

	cout << "Open Hubs for the competitor are==========================" << endl;
	for (int k = 0; k < N; k++)
	{
		if (compY_val[k] > 0)
		{
			cout << k + 1 << endl;
		}
	}

	//====================================Entrant's problem================================================================//


	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			R_val[i][j] = 0;
			Y_val[i] = 0;
		}
	}

	IloNum UB = IloInfinity;
	IloNum best_UB = IloInfinity;
	IloNum LB = 0;
	IloNum best_LB = 0;
	IloNum last_LB = 0;

	//Begin Iterations==============================================================================

	int iter = 0;

	//IloNumVar theta(env, 0, IloInfinity, ILOFLOAT);
	IloExpr Obj(env); // Creates an expression with the name Obj (Objective)
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			Obj += h[i][j] * Z[i][j];
		}
	}


	//Constraint 1:
	IloExpr SumYp(env);

	for (int k = 0; k<N; k++)

	{

		SumYp += Y[k];

	}

	model.add(SumYp == p);

	SumYp.end();

	//Constraint 2:
	for (int i = 0; i<N; i++)

	{

		for (int j = 0; j<N; j++)

		{
			IloExpr SumX(env);
			for (int m = 0; m < N; m++)

			{

				for (int k = 0; k < N; k++)

				{

					SumX += X[i][j][k][m];


				}

			}
			model.add(SumX == 1);
			SumX.end();

		}

	}




	//Constraint 3: for Facet Inducing

	for (int i = 0; i<N; i++)

	{

		for (int j = 0; j<N; j++)

		{

			for (int k = 0; k<N; k++)

			{

				IloExpr SumX(env);

				for (int m = 0; m < N; m++)

				{
					SumX += X[i][j][k][m];

					if (m != k)
					{
						SumX += X[i][j][m][k];
					}
				}

				model.add(SumX <= Y[k]);
				SumX.end();

			}

		}

	}


	//Constraint 5:

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			IloExpr Sum_o(env);
			for (int m = 0; m < N; m++)
			{
				for (int n = 0; n < N; n++)
				{
					Sum_o += X[i][j][m][n] * o[i][j][m][n];
				}
			}
			model.add(R[i][j] == Sum_o);
			Sum_o.end();
		}
	}

	IloNumArray Cuts_approx_init(env, 32, 0.0, 0.0326554, 0.102376, 0.179404, 0.264797, 0.359813, 0.465954, 0.585027, 0.719222, 0.871213,
		1.04429, 1.24255, 1.47111, 1.7365, 2.04706, 2.41367, 2.85069, 3.37736, 4.02001, 4.8154,
		5.81609, 7.09939, 8.78276, 11.0518, 14.2139, 18.8083, 25.854, 37.4721, 58.7112, 104.244,
		233.952, 988.484);

	IloInt H = Cuts_approx_init.getSize();

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			IloNum Sum_s = 0;
			for (int m = 0; m < N; m++)
			{
				for (int n = 0; n < N; n++)
				{
					Sum_s += s[i][j][m][n] * comp_X_val[i][j][m][n];
				}
			}

			for (IloInt h = H - 1; h < H; h++)
			{
				model.add(Z[i][j] <= (pow(Cuts_approx_init[h], 2) + R[i][j] * Sum_s) / (pow((Cuts_approx_init[h] + Sum_s), 2)));
			}
		}
	}

	start = clock();
	// model.add is a function to add constraints and objectives in the CPLEX environment
	model.add(IloMaximize(env, Obj)); // IloMinimize is used for minimization problems
	Obj.end(); // Clear out the memory allocated for the expression 
	// or is ||
	while (((1 - best_LB / best_UB) > .01) && run_SB_tot < TIMLIM)
	{
		start_iter = clock();
		iter += 1;
		cout << endl << "Iteration " << iter << endl;
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				IloNum Sum_s = 0;
				for (int m = 0; m < N; m++)
				{
					for (int n = 0; n < N; n++)
					{
						Sum_s += s[i][j][m][n] * comp_X_val[i][j][m][n];
					}
				}
				model.add(Z[i][j] <= (pow(R_val[i][j], 2) + R[i][j] * Sum_s) / (pow((R_val[i][j] + Sum_s), 2)));
			}
		}
		cplex.setParam(IloCplex::TiLim, TIMLIM);
		// Print results
		IloCplex cplex_(model);
		//cplex.setOut(env.getNullStream()); // if we get an empty stream in return
		//cplex.setWarning(env.getNullStream());

		cplex.solve();//solving the MODEL

		if (cplex.getStatus() == IloAlgorithm::Infeasible) // if the problem is infeasible
		{
			cout << "Problem Infeasible" << endl;
		}
		UB = cplex.getObjValue();

		run_SB = clock() - start_iter;
		// Assign W values to W_val

		//cout << "W----------------------------------" << endl;
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				for (int m = 0; m < N; m++)
				{
					for (int n = 0; n < N; n++)
					{
						X_val[i][j][m][n] = cplex.getValue(X[i][j][m][n]);
					}
				}
			}
		}

			for (int k = 0; k < N; k++)
			{
				Y_val[k] = cplex.getValue(Y[k]);
			}

		IloNum Z_LB = 0;

		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				IloNum Sum_s = 0;
				IloNum Sum_o = 0;
				for (int m = 0; m < N; m++)
				{
					for (int n = 0; n < N; n++)
					{
						Sum_s += s[i][j][m][n] * comp_X_val[i][j][m][n];
						Sum_o += X_val[i][j][m][n] * o[i][j][m][n];
					}
				}

				R_val[i][j] = Sum_o;
				Z_LB += h[i][j] * (R_val[i][j] / (Sum_s + R_val[i][j]));
			}
		}

		last_LB = Z_LB;
		best_LB = IloMax(best_LB, Z_LB);
		best_UB = IloMin(best_UB, UB);

		if ((best_LB - last_LB)<0.01)
		{
			for (int k = 0; k < N; k++)
			{
				best_Y[k] = Y_val[k];
			}
		}

		cout << "The Upper bound on the SACOHLP is " << best_UB << endl;
		cout << "The Lower bound on the SACOHLP is " << best_LB << endl;
		cout << "Run time solving the sub-problem is " << run_SB / (double)CLOCKS_PER_SEC << endl;
		run_SB_tot = run_SB_tot + run_SB;

	}//end while

	
	cout << "==========================Original problem==================================" << endl;
	cout << "The Lower bound to the SACOHLP is " << best_LB << endl;
	cout << "The Upper bound to the SACOHLP is " << best_UB << endl;

	//system("Pause");
	//End of setting and solving SUB2--------------------------------------------------------------------------


	cout << "Optimal hub locations are==========================" << endl;
	for (int k = 0; k < N; k++)
	{
		if (best_Y[k] > eps)
		{
			cout << k + 1 << endl;
		}
	}
	cout << "Total Run time \t" << (clock() - start) / (double)CLOCKS_PER_SEC << endl;
	cout << "Total Run time for sub-problems \t" << run_SB_tot / (double)CLOCKS_PER_SEC << endl;

}

void LR_CPA(int p, double al) //Refer to Section 5, Lagrangian Relaxation with Cutting Plane Algorithm (LR-CPA)
{
	IloInt A;//Number of nodes and hubs
	IloNum eps;//later assigned as eps = cplex.getParam(IloCplex::EpInt); EpInt is the tolerance gap for integer variables
	IloNum be, ga, de;
	A = 1000;
	IloNum c_tranship = 1;
	IloNum c_distribute = 1;
	be = de = 1;
	ga = .75;

	IloModel model(env);
	IloCplex cplex(model);
	Num4DMatrix cost_xy(env, N);
	Num4DMatrix time_xy(env, N);
	Num4DMatrix s(env, N);
	Num4DMatrix o(env, N);

	////////////////////////DEFINE VARIABLES///////////////////////

	for (int i = 0; i < N; i++)
	{
		cost_xy[i] = Num3DMatrix(env, N);
		time_xy[i] = Num3DMatrix(env, N);
		s[i] = Num3DMatrix(env, N);
		o[i] = Num3DMatrix(env, N);
		for (int j = 0; j < N; j++)
		{
			s[i][j] = Num2DMatrix(env, N);
			o[i][j] = Num2DMatrix(env, N);
			cost_xy[i][j] = Num2DMatrix(env, N);
			time_xy[i][j] = Num2DMatrix(env, N);

			for (int k = 0; k < N; k++)
			{
				s[i][j][k] = IloNumArray(env, N);
				o[i][j][k] = IloNumArray(env, N);
				cost_xy[i][j][k] = IloNumArray(env, N);
				time_xy[i][j][k] = IloNumArray(env, N);
			}
		}
	}


	cout << "Attractiveness index is " << A << endl;
	cout << "cost decay parameter is " << be << endl;
	cout << "time decay parameter is " << de << endl;

	eps = cplex.getParam(IloCplex::EpInt);

	for (int i = 0; i<N; i++)
	{
		for (int j = 0; j<N; j++)
		{
			for (int k = 0; k<N; k++)
			{
				for (int l = 0; l<N; l++)
				{
					if (i != j)
					{
						cost_xy[i][j][k][l] = c_distribute*dist_xy[i][k] + al*dist_xy[k][l] + c_tranship*dist_xy[l][j];
						time_xy[i][j][k][l] = dist_xy[i][k] + dist_xy[k][l] + dist_xy[l][j];
						if (k != l)
						{
							s[i][j][k][l] = A / (ga*pow(time_xy[i][j][k][l], be) + (1 - ga)*pow(cost_xy[i][j][k][l], de));
							o[i][j][k][l] = A / (ga*pow(time_xy[i][j][k][l], be) + (1 - ga)*pow(cost_xy[i][j][k][l], de));
						}
						else if (k == l)
						{
							s[i][j][k][l] = A*1.25 / (ga*pow(time_xy[i][j][k][l], be) + (1 - ga)*pow(cost_xy[i][j][k][l], de));
							o[i][j][k][l] = A*1.25 / (ga*pow(time_xy[i][j][k][l], be) + (1 - ga)*pow(cost_xy[i][j][k][l], de));
						}
					}
					if (i == j)
					{
						cost_xy[i][j][k][l] = 0;
						time_xy[i][j][k][l] = 0;
						s[i][j][k][l] = eps;
						o[i][j][k][l] = eps;
					}
				}
			}
		}
	}

	//Define decision variable
	NumVar4DMatrix X(env, N);//1, if hubs are located at both k and l
	Num4DMatrix X_val(env, N);//1, if hubs are located at both k and l
	IloNumArray Y_val(env, N);
	IloBoolVarArray Y(env, N);
	IloNumArray best_Y(env, N);
	Num2DMatrix Yal(env, N);
	IloBoolVarArray comp_Y(env, N);
	IloNumArray compY_val(env, N);
	NumVar4DMatrix comp_X(env, N);//1, if hubs are located at both k and l
	Num4DMatrix comp_X_val(env, N);//1, if hubs are located at both k and l
	NumVar2DMatrix Z(env, N);
	Num2DMatrix Z_val(env, N);
	NumVar2DMatrix R(env, N);
	Num2DMatrix R_new(env, N);
	Num2DMatrix R_val(env, N);
	Num3DMatrix la(env, N);//lagrange multiplier
	Num3DMatrix best_la(env, N);//best lagrange multiplier

	for (int i = 0; i < N; i++)
	{
		comp_X[i] = NumVar3DMatrix(env, N);
		comp_X_val[i] = Num3DMatrix(env, N);
		for (int j = 0; j < N; j++)
		{
			comp_X[i][j] = NumVar2DMatrix(env, N);
			comp_X_val[i][j] = Num2DMatrix(env, N);
			for (int k = 0; k < N; k++)
			{
				comp_X[i][j][k] = IloNumVarArray(env, N, 0, 1, ILOFLOAT);
				comp_X_val[i][j][k] = IloNumArray(env, N);
			}
		}
	}

	IloModel model_comp(env);
	IloCplex cplex_comp(model_comp);
	IloExpr Obj_comp(env); // Creates an expression with the name Obj (Objective)

	for (int i = 0; i<N; i++)
	{
		for (int j = 0; j<N; j++)
		{
			for (int k = 0; k<N; k++)
			{
				for (int m = 0; m<N; m++)
				{
					Obj_comp += (c_distribute*dist_xy[i][k] + al*dist_xy[k][m] + c_tranship*dist_xy[m][j]) *comp_X[i][j][k][m];
				}
			}
		}
	}


	model_comp.add(IloMinimize(env, Obj_comp)); // IloMinimize is used for minimization problems
	Obj_comp.end(); // Clear out the memory allocated for the expression 

	//Constraint 1: 

	for (int i = 0; i<N; i++)

	{

		for (int j = 0; j<N; j++)

		{
			IloExpr SumX(env);

			for (int m = 0; m < N; m++)
			{
				for (int k = 0; k < N; k++)

				{
					SumX += comp_X[i][j][k][m];
				}

			}
			model_comp.add(SumX == 1);
			SumX.end();
		}
	}

	//Constraint 2: sum {k in 1..N}(Y[k])=p;

	IloExpr SumZkk(env);

	for (int k = 0; k<N; k++)
	{

		SumZkk += comp_Y[k];

	}

	model_comp.add(SumZkk == p);

	SumZkk.end();



	//Constraint 3: for Facet Inducing

	for (int i = 0; i<N; i++)

	{

		for (int j = 0; j<N; j++)

		{

			for (int k = 0; k<N; k++)

			{

				IloExpr SumX(env);

				for (int m = 0; m < N; m++)

				{
					SumX += comp_X[i][j][k][m];

					if (m != k)
					{
						SumX += comp_X[i][j][m][k];
					}
				}

				model_comp.add(SumX <= comp_Y[k]);
				SumX.end();

			}

		}

	}

	//==============================================================================

	//Optimize

	cplex_comp.setOut(env.getNullStream()); //This is to supress the output of Branch & Bound Tree on screen

	cplex_comp.setWarning(env.getNullStream()); //This is to supress warning messages on screen

	cplex_comp.solve();//solving the MODEL

	if (cplex_comp.getStatus() == IloAlgorithm::Infeasible) // if the problem is infeasible

	{

		cout << "Problem is Infeasible" << endl;

	}

	for (int k = 0; k<N; k++)

	{

		if (cplex_comp.getValue(comp_Y[k]) > 1 - eps)

		{

			compY_val[k] = 1;

		}

		else

		{

			compY_val[k] = 0;

		}

	}

	for (int i = 0; i<N; i++)

	{

		for (int j = 0; j<N; j++)

		{

			for (int k = 0; k<N; k++)

			{

				for (int m = 0; m<N; m++)

				{

					if (cplex_comp.getValue(comp_X[i][j][k][m]) > 1 - eps)

					{

						comp_X_val[i][j][k][m] = 1;

					}

					else

					{

						comp_X_val[i][j][k][m] = 0;

					}

				}

			}

		}

	}

	cout << "Open Hubs for the competitor are==========================" << endl;
	for (int k = 0; k < N; k++)
	{
		if (compY_val[k] > 0)
		{
			cout << k + 1 << endl;
		}
	}

	//====================================Entrant's problem================================================================//

	clock_t start, run_SB, run_iter = 0, run_SB_tot = 0, tot_time;
	start = clock();
	tot_time = clock();
	int Compareval(const void *, const void *);
	XVAL *alpha;
	XVAL *sorted;//for storing the values of x[i][j] variables for sorting by qsort library function
	alpha = (XVAL *)calloc(N*N*N, sizeof(XVAL));
	sorted = (XVAL *)calloc(N, sizeof(XVAL));
	//Assigning value to variables

	for (int i = 0; i < N*N*N; i++)
	{
		alpha[i].index = i;
		alpha[i].value = 0;
	}

	for (int i = 0; i < N; i++)
	{
		la[i] = Num2DMatrix(env, N);
		best_la[i] = Num2DMatrix(env, N);
		for (int j = 0; j < N; j++)
		{
			la[i][j] = IloNumArray(env, N, 0, IloInfinity, ILOFLOAT);
			best_la[i][j] = IloNumArray(env, N, 0, IloInfinity, ILOFLOAT);
		}
	}

	for (int i = 0; i < N; i++)
	{
		X[i] = NumVar3DMatrix(env, N);
		X_val[i] = Num3DMatrix(env, N);
		Z[i] = IloNumVarArray(env, N, 0, 1, ILOFLOAT);
		Z_val[i] = IloNumArray(env, N);
		Yal[i] = IloNumArray(env, N);
		R[i] = IloNumVarArray(env, N, 0, IloInfinity, ILOFLOAT);
		R_val[i] = IloNumArray(env, N);
		R_new[i] = IloNumArray(env, N);
		for (int j = 0; j < N; j++)
		{
			X[i][j] = NumVar2DMatrix(env, N);
			X_val[i][j] = Num2DMatrix(env, N);
			for (int k = 0; k < N; k++)
			{
				X[i][j][k] = IloNumVarArray(env, N, 0, 1, ILOFLOAT);
				X_val[i][j][k] = IloNumArray(env, N);
			}
		}
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			R_val[i][j] = 0;
		}
	}

	IloNum UB = IloInfinity;
	IloNum best_UB = IloInfinity;
	IloNum LB = 0;
	IloNum best_LB = 0;

	//quotient = dividend / divisor;
	//remainder = dividend % divisor;
	//eps = cplex.getParam(IloCplex::EpInt);
	//cplex.setParam(IloCplex::EpAGap, 0.01);
	//cplex.setParam(IloCplex::TiLim, 14400);
	//clock_t t1, t2;
	//t1 = clock();
	//cplex.setParam(IloCplex::MIPEmphasis, 2);
	//cplex.solve();//solving the MODEL
	//t2 = clock();

	//Begin Iterations==============================================================================

	IloNum t = 0;
	IloNum del = 2;
	int iter = 0;


	while ((1 - best_LB / best_UB) > .01)
	{
		run_iter = clock();
		iter += 1;
		cout << endl << "Iteration " << iter << endl;

		for (int i = 0; i < N; i++)
		{
			for (int k = 0; k < N; k++)
			{
				Y_val[k] = 0;
				Yal[i][k] = 0;
			}

		}
		
		//=============Set-up and solve SUB1 without solver============================================================
		
		int count = 0;
		IloNum Sumsort;
		for (int k = 0; k < N; k++)
		{
			sorted[count].index = k;
			sorted[count].value = 0;
			
			for (int i = 0; i < N; i++)
			{
				Sumsort = 0;
				for (int j = 0; j < N; j++)
				{
					Sumsort += alpha[N*N*i + N*j + k].value ;
				}
				sorted[count].value += Sumsort;
			}
			count++;
		}
		
		qsort((XVAL *)sorted, count, sizeof(sorted[0]), Comparevalue);

	https://www.tutorialspoint.com/c_standard_library/c_function_qsort.htm


		for (int m = 0; m < p; m++)
		{
			Y_val[(sorted[m].index)] = 1;
		}

		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				for (int k = 0; k < N; k++)
				{
					la[i][j][k] = alpha[N*N*i + N*j + k].value;
					best_la[i][j][k] = alpha[N*N*i + N*j + k].value;
				}
			}
		}
		
		cout << "Open Hubs are==========================" << endl;
		for (int k = 0; k < N; k++)
		{
			if (Y_val[k] > 0)
			{
				cout << k + 1 << endl;
			}
		}

		// Generating a feasible solution

		IloNumArray Clos(env, N);
		for (int i = 0; i < N; i++)
		{
			Clos[i] = IloInfinity;
			for (int k = 0; k < N; k++)
			{
				if (Y_val[k] == 1)
				{
					Clos[i] = IloMin(Clos[i], dist_xy[i][k]);
				}
				
			}

			for (int k = 0; k < N; k++)
			{
				if (Y_val[k] == 1)
				{
					if (Clos[i] == dist_xy[i][k])
					{
						Yal[i][k] = 1;
					}
				}
			}
				
		}
		
		//cout << "=================Operational routes are==================" << endl;
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				for (int k = 0; k < N; k++)
				{
					for (int l = 0; l < N; l++)
					{
						if (Yal[i][k] * Yal[j][l] > 0)
						{
							//cout << i + 1 << "--->" << k + 1 << "--->" << l + 1 << "--->" << j + 1 << endl;
						}
					}
				}
			}
		}

		IloNum SumLB = 0;
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				IloNum Suma = 0;
				IloNum Sumx = 0;
				for (int k = 0; k < N; k++)
				{
					for (int l = 0; l < N; l++)
					{
						Suma += s[i][j][k][l] * comp_X_val[i][j][k][l];
						Sumx += o[i][j][k][l] *  Yal[i][k]*Yal[j][l];
					}
				}
				SumLB += h[i][j] * (Sumx / (Suma + Sumx));
			}
		}

		IloNum last_LB = LB;
		LB = IloMax(LB, SumLB);


		int iteration = 0;
		IloNum UB_SUB2 = 1e12;
		IloNum LB_SUB2 = -1e12;

		//=============================================================================================
		//Set up and solve SUB-2-----------------------------------------------------------------------
		//=============================================================================================

		IloModel model_SUB2(env);
		//IloNumVar theta(env, 0, IloInfinity, ILOFLOAT);
		IloExpr Obj_SUB2(env); // Creates an expression with the name Obj (Objective)

		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				Obj_SUB2 += h[i][j] * Z[i][j];
				for (int k = 0; k < N; k++)
				{
					IloExpr Suml1(env);
					IloExpr Suml2(env);
					for (int l = 0; l < N; l++)
					{
						Suml1 += X[i][j][k][l];

						if (k != l)
						{
							Suml2 += X[i][j][l][k];
						}
					}
					Obj_SUB2 -= la[i][j][k] * (Suml1+Suml2);
				}
			}
		}
	

		//Constraint 1:
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				IloExpr Sum(env);
				for (int k = 0; k < N; k++)
				{
					for (int l = 0; l < N; l++)
					{
						Sum += X[i][j][k][l];
					}
				}
				model_SUB2.add(Sum == 1);
				Sum.end();
			}
		}

		//Constraint 5:

		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				IloExpr Sum_o(env);
				for (int m = 0; m < N; m++)
				{
					for (int n = 0; n < N; n++)
					{
						Sum_o += X[i][j][m][n] * o[i][j][m][n];
					}
				}
				model_SUB2.add(R[i][j] == Sum_o);
				Sum_o.end();
			}
		}

		// model.add is a function to add constraints and objectives in the CPLEX environment
		model_SUB2.add(IloMaximize(env, Obj_SUB2)); // IloMinimize is used for minimization problems
		Obj_SUB2.end(); // Clear out the memory allocated for the expression 
		IloNumArray Cuts_approx_init(env, 32, 0.0, 0.0326554, 0.102376, 0.179404, 0.264797, 0.359813, 0.465954, 0.585027, 0.719222, 0.871213,
			1.04429, 1.24255, 1.47111, 1.7365, 2.04706, 2.41367, 2.85069, 3.37736, 4.02001, 4.8154,
			5.81609, 7.09939, 8.78276, 11.0518, 14.2139, 18.8083, 25.854, 37.4721, 58.7112, 104.244,
			233.952, 988.484);

		IloInt H = Cuts_approx_init.getSize();
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				IloNum Sum_s = 0;
				for (int m = 0; m < N; m++)
				{
					for (int n = 0; n < N; n++)
					{
						Sum_s += s[i][j][m][n] * comp_X_val[i][j][m][n];
					}
				}
				for (IloInt h = H - 1; h < H; h++)
				{
					model.add(Z[i][j] <= (pow(Cuts_approx_init[h], 2) + R[i][j] * Sum_s) / (pow((Cuts_approx_init[h] + Sum_s), 2)));
				}
			}
		}

		run_SB = clock();
		while (1 - (LB_SUB2 / UB_SUB2) > .01)
		{
			run_SB = clock();
			iteration += 1;
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					IloNum Sum_s = 0;
					for (int m = 0; m < N; m++)
					{
						for (int n = 0; n < N; n++)
						{
							Sum_s += s[i][j][m][n]*comp_X_val[i][j][m][n];
						}
					}
					model_SUB2.add(Z[i][j] <= (pow(R_val[i][j], 2) + R[i][j] * Sum_s) / (pow((R_val[i][j] + Sum_s), 2)));
				}
			}

			// Print results
			IloCplex cplex_SUB2(model_SUB2);
			cplex_SUB2.setOut(env.getNullStream()); // if we get an empty stream in return
			cplex_SUB2.setWarning(env.getNullStream());

			cplex_SUB2.solve();//solving the MODEL

			if (cplex_SUB2.getStatus() == IloAlgorithm::Infeasible) // if the problem is infeasible
			{
				cout << "SUB2 Infeasible" << endl;
			}
			UB_SUB2 = cplex_SUB2.getObjValue();
			//cout << "The upper bound to the sub-problem at iteration" << " " << iteration << UB_SUB2 <<endl;

			run_SB = clock() - run_SB;
			// Assign W values to W_val

			//cout << "W----------------------------------" << endl;
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					Z_val[i][j] = cplex_SUB2.getValue(Z[i][j]);
					for (int m = 0; m < N; m++)
					{
						for (int n = 0; n < N; n++)
						{
							X_val[i][j][m][n] = cplex_SUB2.getValue(X[i][j][m][n]);

						}
					}
				}
			}


			IloNum Z_LB = 0;

			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					IloNum Sum_s = 0;
					IloNum Sum_o = 0;
					for (int m = 0; m < N; m++)
					{
						for (int n = 0; n < N; n++)
						{
							Sum_s += s[i][j][m][n]*comp_X_val[i][j][m][n];
							Sum_o += o[i][j][m][n]*X_val[i][j][m][n];
						}
					}
					R_val[i][j] = Sum_o;
					Z_LB += h[i][j] * (R_val[i][j] / (Sum_s + R_val[i][j]));
				}
			}
			

			
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					for (int k = 0; k < N; k++)
					{
						IloNum Suml1 = 0;
						IloNum Suml2 = 0;
						for (int l = 0; l < N; l++)
						{
							Suml1 += X_val[i][j][k][l];

							if (k != l)
							{
								Suml2 += X_val[i][j][l][k];
							}
						}
						Z_LB = Z_LB - (la[i][j][k])*(Suml1 + Suml2);
					}
				}
			}

			LB_SUB2 = IloMax(LB_SUB2, Z_LB);
			run_SB_tot += run_SB;
			cout << "The Upper bound on the sub problem is " << UB_SUB2 << endl;
			cout << "The Lower bound on the sub problem is " << LB_SUB2 << endl;
		}//end while
	
		cout << "The Upper bound on the sub problem is " << UB_SUB2 << endl;
		cout << "The Lower bound on the sub problem is " << LB_SUB2 << endl;
		cout << "Run time solving the sub-problem is " << run_SB / (double)CLOCKS_PER_SEC << endl;
		run_SB_tot += run_SB;

		IloNum theta_SUB1 = 0;
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				for (int k = 0; k < N; k++)
				{
					theta_SUB1 += (la[i][j][k]) * Y_val[k];
				}
			}
		}

		IloNum last_UB = UB;
		UB = UB_SUB2 + theta_SUB1;
		best_UB = IloMin(best_UB, UB);

		IloNum Sums = 0;
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				for (int k = 0; k < N; k++)
				{
					IloNum Sum1 = 0;
					IloNum Sum2 = 0;
					for (int l = 0; l < N; l++)
					{
						Sum1 += X_val[i][j][k][l];
						if (k != l)
						{
							Sum2 += X_val[i][j][l][k];
						}
					}
					Sums += pow((Y_val[k] - Sum1 - Sum2), 2);
				}
			}
		}

		t = del*(best_UB - best_LB) / Sums;

		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				for (int k = 0; k < N; k++)
				{
					IloNum Sum1 = 0;
					IloNum Sum2 = 0;
					for (int l = 0; l < N; l++)
					{
						Sum1 += X_val[i][j][k][l];
						if (k != l)
						{
							Sum2 += X_val[i][j][l][k];
						}
					}
				
					la[i][j][k] = la[i][j][k] - t*(Y_val[k] - Sum1 - Sum2);
						
				}
			}
		}

		
		
		if ((LB - last_LB)<0.01)
		{
			for (int k = 0; k < N; k++)
			{
				best_Y[k] = Y_val[k];
			}
		}

		if (UB < best_UB)
		{
			best_la = la;
		}

		if ((last_UB - best_UB)<0.01)
		{
			count = count + 1;
		}

		else
		{
			count = 0;
		}
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				for (int k = 0; k < N; k++)
				{
					alpha[N*N*i + N*j + k].value = la[i][j][k];
					if (count == 50)
					{
						alpha[N*N*i + N*j + k].value = IloMax(0,best_la[i][j][k]);
						del = del / 2;
						cout << "------------------------------Delta Reduced-------------------------------------";
						count = 0;
					}
				}
			}
		}

		
		best_LB = IloMax(LB, best_LB);
		//cout << "current delta: " << del<<endl;
		if (del < 0.1)
		{
			del = 1;
		}
		cout << "The Lower bound to the COHLP at this iteration is " << best_LB << endl;
		cout << "The Upper bound to the COHLP at this iteration is " << best_UB << endl;
		cout << "Run time for iteration \t" << iter << " is " << (clock() - run_iter) / (double)CLOCKS_PER_SEC << endl;

		//End of setting and solving SUB2--------------------------------------------------------------------------
		
		tot_time = clock() - start;
	}//end while

	cout << "Optimal hub locations are==========================" << endl;
	for (int k = 0; k < N; k++)
	{
		if (best_Y[k] > eps)
		{
			cout << k + 1 << endl;
		}
	}

	IloNum MS = 0;
	IloNum Flow = 0;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			Flow = Flow + h[i][j];
		}

	}
	MS = best_LB / Flow;

	cout << "==========================Original problem==================================" << endl;
	cout << "The Lower bound to the COHLP is " << best_LB << endl;
	cout << "The Upper bound to the COHLP is " << best_UB << endl;

	cout << "The market share captured by the entrant is " << MS << endl;
	cout << "Optimality gap for COHLP is " << (1 - best_LB / best_UB) * 100 << " %" << endl;
	cout << "Total Run time \t" << tot_time / (double)CLOCKS_PER_SEC << endl;
	cout << "Total Run time for sub-problems \t" << run_SB_tot / (double)CLOCKS_PER_SEC << endl;

}


int p, alpha;
int main(int argc, char **argv)
{
	try
	{
	const char* data_filename = "Deafult_FileName.richa";
	if (argc > 1)
	{
		data_filename = argv[1];
		p = atoi(argv[2]);
		alpha = atof(argv[3]);
	}
	fstream datafile;
	datafile.open(data_filename, ios::in);

	if (!datafile)
	{
		cerr << "ERROR: could not open file " << data_filename << " for reading" << endl;
		cerr << "usage:   " << argv[0] << " <datafile>" << endl;
		throw(-1);
	}
	datafile >> N >> dist_xy >> h;
	datafile.close();
	std::ofstream outfile;
	outfile.open("Experiment_Results.txt", std::ios_base::app);

	cout << "=============================The SOCP results are" << endl;
	SOCP(p, alpha);

	cout << "=============================The CPA results of Hamacher are" << endl;
	CPA(p, alpha);

	

	cout << "=============================The LR-CPA results of Hamacher are" << endl;
	LR_CPA(p, alpha);
	
	
}

	catch (IloException& ex)
	{
		cerr << "Error: " << ex << endl;
	}
	catch (...)
	{
		cerr << "Error" << endl;
	}

	return 0;
}