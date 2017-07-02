#include <iostream>
#include <deque>
#include <cstdlib>
#include <ctime>
#include <cmath>
#define POP_SIZE 100
//TEST SYSTEM 1, CASE 1


//issues - actual values of P_smin

//Constants
//[main]
double T=24; //time interval
double iterators = 30; //iterations?
double N_h = 4; //number of hydro plants
double N_s = 1; //number of thermal plants
int N_p = 100; //total number of population
int P_e = 10; //total number of elite individuals in N_p
int P_c = 40; //population from crossover
int P_m = 50; //population number of mutated individuals
double e_vcourse = 1; //for constraints handling
double e_vfine = .01;//for constraints handling
double e_pcourse = 2;//
double e_pfine = .01;
double iteration_course = 20;
double iteration_fine = 10;
//[37]
double P_smin = 500;
double P_smax = 2500;
double Q_hmin[4] = {5,6,10,13};
double Q_hmax[4] = {15,15,30,25};
double V_hmin[4] = {80,60,100,170};
double V_hmax[4] = {150,120,240,160};
double V_hbegin[4] = {100,80,170,120};
double V_hend[4] = {120,70,170,140};
//[main]
double a_si = 5000;
double b_si = 19.2;
double c_si = .002;
double d_si = 700;
double e_si = .085;
double UR[2] = {0,0};
double DR[2] = {0,0};
//[37]
double R_u[4] = {0,0,2,1};
double P_D[24] = {1370, 1390, 1360, 1290, 1290, 1410, 1650, 2000, 2240, 2320, 2230, 2310, 2230, 2200, 2130, 2070, 2130, 2140, 2240, 2280, 2240, 2120, 1850, 1590};
double time_delay[4] = {2,3,4,0};
//[main]
double gmax_rcga = 1000;
double gmax_pso = 400;
//[39]
double b_mutation = 1;
double n_c = 1;
// P_s = zeros()
// P_h =zeros()

// P_L = zeros()
// F = null
//[37]
double c1[4] = {-0.0042, -0.0040, -0.0016, -0.0030};
double c2[4] = {-0.42, -0.30, -0.30, -0.31};
double c3[4] = {0.030, .015, 0.014, 0.027};
double c4[4] = {0.90, 1.14, 0.55, 1.44};
double c5[4] = {10, 9.5, 5.5, 14};
double c6[4] = {-50, -70, -40, -90};

double inflow[4][24];

int upstreamIndex[4][4];







using namespace std;

deque<deque<deque<double> > > population;// the last element in 3d deque is the fitness value
deque<double> population_fitness;
deque<deque<deque<double> > > volume;
deque<deque<deque<double> > > population_elite;
deque<deque<deque<double> > > population_crossover;
deque<deque<deque<double> > > population_mutation;
deque<deque<deque<double> > > velocity;
deque<deque<deque<double> > > population_pso;
deque<double> population_pso_fitness;
deque<deque<deque<double> > > pbest;//equal to its current fitness function
deque<double> pbest_fitness;
deque<deque<double> > gbest;
double gbest_fitness;

double best_fitness;


//SUPER UNIVERSAL VARIABLES
deque<double> aveBestFitness;
deque<double> aveWorstFitness;
deque<deque<double> > aveBestChromosome;
deque<deque<double> > aveWorstChromosome;
deque<deque<double> > aveBestPh;
deque<deque<double> > aveBestPs;
deque<deque<double> > aveBestQh;
double aveTime;
//FUNCTIONS
void inflowF(){//[37]
	double huehue1[24] = {10, 9, 8, 7, 6, 7, 8, 9, 10, 11, 12, 10, 11, 12, 11, 10, 9, 8, 7, 6, 7, 8, 9, 10};
	double huehue2[24] = {8, 8, 9, 9, 8, 7, 6, 7, 8, 9, 9, 8, 8, 9, 9, 8, 7, 6, 7, 8, 9, 9, 8, 8};
	double huehue3[24] = {8.1, 8.2, 4, 2, 3, 4, 3, 2, 1, 1, 1, 2, 4, 3, 3, 2, 2, 2, 1, 1, 2, 2, 1, 0};
	double huehue4[24] = {2.8, 2.4, 1.6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

	for(int i=0;i<24;i++){
		inflow[0][i] = huehue1[i];
	}
	for(int i=0;i<24;i++){
		inflow[1][i] = huehue2[i];
	}
	for(int i=0;i<24;i++){
		inflow[2][i] = huehue3[i];
	}
	for(int i=0;i<24;i++){
		inflow[3][i] = huehue4[i];
	}
	
	
}

void initUpstream(){//initialize what plant is upstream of the jth hydro plant, if -1, no plant upstream
	for(int i=0;i<N_h;i++){
		for(int j=0;j<N_h;j++){
			upstreamIndex[i][j] = -1;
		}
	}
	upstreamIndex[2][0] = 0;
	upstreamIndex[2][1] = 1;
	upstreamIndex[3][0] = 2;
}

double sum_upstream(int pop,int t,int j){
	double temp=0;
	for(int i=0;i<R_u[j];i++){
		if((t-time_delay[i])<0){
			continue;
		}
		else{
			if(isinf(population[pop][t-time_delay[i]][j]) || isinf(population[pop][t-time_delay[i]][j])!=isinf(population[pop][t-time_delay[i]][j])){
				population[pop][t-time_delay[i]][j] = Q_hmax[j];
				continue;
			}
			else{////////////////////////////////////////////
					temp = temp + population[pop][t-time_delay[i]][upstreamIndex[j][i]];
				
				//temp = temp + population[pop][t-time_delay[i]][j];
			}
			
			//cout<<"sumup: "<<population[pop][t-time_delay[i]][j]<<endl;
		}
	}
	//cout<<"Individual: "<<pop<<endl;
	//cout<<t<<" "<<j<<" "<<temp<<endl;
	return temp;
}

void check_discharge(int index,int j){//check every Q_h for all t in jth hydro plant
	for(int t=0;t<T;t++){
		if(population[index][t][j]<Q_hmin[j]){
			population[index][t][j] = Q_hmin[j];
			if(population[index][t][j]!=population[index][t][j] || isinf(population[index][t][j])){
				population[index][t][j] = Q_hmax[j];
				// cout<<"check_discharge min"<<endl;
				// exit(EXIT_FAILURE);
			}
		}
		if(population[index][t][j]>Q_hmax[j]){
			population[index][t][j] = Q_hmax[j];
			if(population[index][t][j]!=population[index][t][j] || isinf(population[index][t][j])){
				population[index][t][j] = Q_hmax[j];
				// cout<<"check_discharge max"<<endl;
				// exit(EXIT_FAILURE);
			}
		}
	}
}

void check_discharge1(int index,int j, int t){//check every Q_h for all t in jth hydro plant
		if(population[index][t][j]<Q_hmin[j]){
			// cout<<"check_discharge1 less"<<endl;
			// cout<<population[index][t][j]<<endl;
			population[index][t][j] = Q_hmin[j];
			if(population[index][t][j]!=population[index][t][j] || isinf(population[index][t][j])){
				// cout<<"check_discharge1 min"<<endl;
				population[index][t][j] = Q_hmax[j];
				// exit(EXIT_FAILURE);
			}
		}
		if(population[index][t][j]>Q_hmax[j]){
			population[index][t][j] = Q_hmax[j];
			if(population[index][t][j]!=population[index][t][j] || isinf(population[index][t][j])){
				population[index][t][j] = Q_hmax[j];
				// cout<<"check_discharge1 max"<<endl;
				// exit(EXIT_FAILURE);
			}
		}
}

void check_thermal1(int index, int t, int j){
	if(population[index][t][j] <=P_smin){
		population[index][t][j] = P_smin;
		if(population[index][t][j]!=population[index][t][j]){
				//cout<<"check_thermal1 min"<<endl;
				//exit(EXIT_FAILURE);
			}
	}
	if(population[index][t][j] >= P_smax){
		population[index][t][j] = P_smax;
		if(population[index][t][j]!=population[index][t][j]){
				//cout<<"check_thermal1 max"<<endl;
				//exit(EXIT_FAILURE);
			}
	}
}

deque<double> Q_hcF(int index,int j){
	deque<double> temp;
	for(int t=0;t<T;t++){
		if(t==0){
			temp.push_back(population[index][t][j] - Q_hmin[j]);
			//cout<<population[index][t][j] - Q_hmin[j]<<endl;
		}
		else{
			temp.push_back((population[index][t][j]-population[index][t-1][j]));
			//cout<<(population[index][t][j]-population[index][t-1][j])<<endl;
			if(population[index][t][j]!=population[index][t][j] || population[index][t-1][j]!=population[index][t-1][j] || isinf(population[index][t][j]) || isinf(population[index][t-1][j])){
				population[index][t][j] = Q_hmax[j];
				population[index][t-1][j] = Q_hmax[j];
				// cout<<"Q_hcF:["<<t<<"] ["<<j<<"]"<<population[index][t][j]<<" "<<population[index][t-1][j]<<endl;
				// cout<<"FAILURE Q_hcF"<<endl;
				// exit(EXIT_FAILURE);
			}

		}
		
	}
	return temp;
}

int max_adjustableQ(deque<double> input){
	int itr=0;
	double temp=99999999999999999;
	for(int i=1;i<T;i++){
		if(temp<input[i]){
			temp = input[i];
			itr=i;
		}
	}
	return itr;
}

int max_adjustableP(int index, int t){
	/* double temp = 99999999999999999;
	int adIndex=0;
	for(int itr1=N_h;itr1<N_h+N_s;itr1++){
		if(population[index][t][itr1]<temp){
			adIndex = itr1;
		}
	}
	return adIndex; */
	return 0;
	/* int itr=0;
	double temp=0;
	for(int ns=N_h;ns<N_h+N_s;ns++){
		if(temp<input[ns]){
			temp = input[ns];
			itr=ns;
		}
	}
	return itr; */
}

double P_ssum(int index, int t){
	double temp=0;
	for(int a = N_h; a<N_h+N_s;a++){
		if(population[index][t][a]<0){
			continue;
		}
		temp = temp + population[index][t][a];
		
	}
	/* if(temp!=temp){
		cout<<"P_ssum"<<endl;
		//exit(EXIT_FAILURE);
	} */
	if(temp<0){
		temp =0;
	}
	return temp;
}

double P_hsum(int index, int t){//work on the equation
	double temp=0;
	for(int a = 0; a<N_h;a++){
		double temp1 = c1[a]*pow(volume[index][t][a],2) + c2[a]*pow(population[index][t][a],2) + c3[a]*volume[index][t][a]*population[index][t][a] + c4[a]*volume[index][t][a] + c5[a]*population[index][t][a] + c6[a];
		if(temp1<0){
			continue;
		}
		temp = temp + temp1;
	}
	//cout<<index<<" "<<t<<endl;
	//cout<<"Hydropower: "<<temp<<endl;
	//cout<<endl;
	/* if(temp!=temp){
		cout<<"P_ssum"<<endl;
		//exit(EXIT_FAILURE);
	} */
	if(temp<0){
		temp = 0;
	}
	return temp;
}

void update_volume(){
	for(int i=0;i<POP_SIZE;i++){
		deque<deque<double> > V_h;//V_h(t,j)
		V_h.clear();
		for(int t=0;t<T;t++){
			deque<double> temp1;
			for(int j=0;j<N_h;j++){
				double temp=0;//V_h(j)
				if((t-1)<0){
					temp = V_hbegin[j];// + inflow[j][t] - population[i][t][j] + sum_upstream(i,t,j);//spillage not calculated
					temp1.push_back(temp);
					/* if(temp!=temp){
						//cout<<"FAILURE update volume1"<<endl;
						exit(EXIT_FAILURE);
					} */
				}
				else{
					temp = V_h[t-1][j] + inflow[j][t] - population[i][t][j] + sum_upstream(i,t,j);//spillage not calculated
					temp1.push_back(temp);
					/* if(temp!=temp){
						//cout<<"FAILURE update volume2"<<endl;
						exit(EXIT_FAILURE);
					} */
				}
			}
			V_h.push_back(temp1);
			temp1.clear();
		}
		volume[i].swap(V_h);
		V_h.clear();
	}
	
}

void update_volume(int index){
	//cout<<"uv start"<<endl;
	//cout<<"index: "<<index<<endl;
	//cout<<"volume size: "<<volume.size()<<endl;
	//cout<<"volume[i] size: "<<volume[index].size()<<endl;
	//cout<<"pop size: "<<population[index].size()<<endl; 
	deque<deque<double> > V_h;//V_h(t,j)
		V_h.clear();
		for(int t=0;t<T;t++){
			deque<double> temp1;
			for(int j=0;j<N_h;j++){
				double temp;//V_h(j)
				if(t==0){
					temp = V_hbegin[j];// + inflow[j][t] - population[index][t][j] + sum_upstream(index,t,j);//spillage not calculated
					temp1.push_back(temp);
					if(isinf(temp) || temp!=temp){
						cout<<"FAILURE update volume1 i"<<endl;
						exit(EXIT_FAILURE);
					}
					//cout<<"check 1"<<endl;
				}
				else{
					temp = V_h[t-1][j] + inflow[j][t] - population[index][t][j] + sum_upstream(index,t,j);//spillage not calculated
					
					temp1.push_back(temp);
					if(isinf(temp)||temp!=temp){
						cout<<"temp "<<temp<<" inflow "<<inflow[j][t]<<" population "<<population[index][t][j]<<" sum_upst "<<sum_upstream(index,t,j)<<endl;
						cout<<"FAILURE update volume2 i"<<endl;
						exit(EXIT_FAILURE);
					}
					//cout<<"check 2"<<endl;
				}
				//cout<<"j: "<<j<<endl;
			}
			V_h.push_back(temp1);
			temp1.clear();
			//cout<<"t: "<<t<<endl;
		}
		if(volume[index].size()!=V_h.size()){
			cout<<"volume and V_h !size"<<endl;
			exit(EXIT_FAILURE);
		}
		//cout<<"volume[index].size() "<<volume[index].size()<<endl;
		//cout<<"V_h.size() "<<V_h.size()<<endl;
		for(int t=0;t<T;t++){
			for(int j=0;j<N_h;j++){
				volume[index][t][j] = V_h[t][j];
			}
		}
		//volume[index].swap(V_h);
		V_h.clear();
	//cout<<"uv end"<<endl;
}

void evaluation(){//calculate fitness value for population
//cout<<"Check!"<<endl;
	population_fitness.clear();
	for(int itr=0;itr<POP_SIZE;itr++){
		double temp=0;
		for(int t=0;t<T;t++){
			for(int i=0;i<N_s;i++){
				//cout<<itr<<" "<<t<<" "<<i<<endl;
				temp = temp + a_si + b_si*population[itr][t][N_h+i] + c_si*pow(population[itr][t][N_h+i],2);// + fabs(d_si*sin(e_si*(P_smin - population[itr][t][N_h+i])));
				//cout<<temp<<" ";
			}
			//cout<<endl;
			
		}

		population_fitness.push_back(temp);//default place for the fitness value
	}
	
	//bubble sort
	for(int itr1=0;itr1<POP_SIZE;itr1++){
		for(int itr2=0;itr2<POP_SIZE;itr2++){
			if(population_fitness[itr1]<population_fitness[itr2]){
				population[itr1].swap(population[itr2]);
				double huehue = population_fitness[itr1];
				population_fitness[itr1] = population_fitness[itr2];
				population_fitness[itr2] = huehue;
				//population_fitness[itr1].swap(population_fitness[itr2]);
			}
		}
	}
}

void elitist(){//get all elite population and best individual
	population_elite.clear();
	if(best_fitness>population_fitness[0]){
		best_fitness = population_fitness[0];
	}
	for(int i=0;i<P_e;i++){
		population_elite.push_back(population[i]);
	}
}

void sbx(){
	double u;
	do{
		u = (double)rand()/rand();
		u-=(int)u;
	}while(isinf(u) || u!=u);
	
	//deque<deque<deque<double> > > temp;
	for(int itr=0;itr<P_c;itr+=2){
		
		for(int t=0;t<T;t++){
		
			for(int i=0;i<N_h;i++){
				//compute for beta,gamma, u
				double x1;
				double x2;
				double min;
				if(population_crossover[itr][t][i]>=population_crossover[itr+1][t][i]){
					x1 = population_crossover[itr][t][i];
					x2 = population_crossover[itr+1][t][i];
				}
				if(population_crossover[itr][t][i]<population_crossover[itr+1][t][i]){
					x1 = population_crossover[itr][t][i];
					x2 = population_crossover[itr+1][t][i];
				}
				if(population_crossover[itr][t][i]!=population_crossover[itr][t][i] || isinf(population_crossover[itr][t][i])){
					population_crossover[itr][t][i] = Q_hmax[i];
					// cout<<population_crossover[itr][t][i]<<endl;
					// cout<<"FAILURE sbx N_h 0"<<endl;
					i--;
					continue;
					//exit(EXIT_FAILURE);
				}
				if(population_crossover[itr+1][t][i]!=population_crossover[itr+1][t][i] || isinf(population_crossover[itr+1][t][i])){
					population_crossover[itr+1][t][i] = Q_hmax[i];
					i--;
					continue;
					// cout<<population_crossover[itr][t][i]<<endl;
					// cout<<"FAILURE sbx N_h 1"<<endl;
					//exit(EXIT_FAILURE);
				}
				if(x1==x2){
					continue;
				}
				else{
					if((x1-Q_hmin[i])<(x2-Q_hmax[i])){
						min = x1-Q_hmin[i];
					}
					else{
						min = x2-Q_hmax[i];
					}
					double beta = 1 + ((double)2/(x2-x1))*min;
						if(x2!=x2 || isinf(x2)){
							i--;
							continue;
							cout<<"beta x2"<<endl;
							exit(EXIT_FAILURE);
						}
						if(x1!=x1 || isinf(x1)){
							i--;
							continue;
							cout<<"beta x1"<<endl;
							exit(EXIT_FAILURE);
						}
					if(beta!=beta || isinf(beta)){
						i--;
						continue;
						cout<<x2<<" "<<x1<<endl;
						cout<<"FAILURE beta!"<<endl;
						exit(EXIT_FAILURE);
					}
					//cout<<-(n_c+1)<<endl;
					double asdf = -(n_c+1);
					//cout<<"asdf "<<asdf<<endl;
					double alpha;
					if(beta<0){
						alpha = 2 - -pow(fabs(beta),asdf);
					}
					else{
						alpha = 2 - pow(fabs(beta),asdf);
					}
					if(alpha!=alpha || isinf(alpha)){
						i--;
						continue;
						cout<<alpha<<" "<<beta<<" "<<asdf<<endl;
						exit(EXIT_FAILURE);
					}
					//double alpha = 2 - pow(beta,asdf);
					/* cout<<"pow beta "<<pow(beta,(-(n_c+1)))<<endl;
					cout<<"n_c "<<n_c<<endl;
					cout<<"beta "<<beta<<endl;
					cout<<"alpha "<<alpha<<endl;
					cout<<"u "<<u<<endl; */
					double gamma = 0;
					if(u<=(1/alpha)){
						if((alpha*u)<0){
							gamma = -pow(fabs((alpha*u)),((double)1/(n_c+1)));
						}
						else{
							gamma = pow((alpha*u),((double)1/(n_c+1)));
						}
					}
					else{
						if((1/(2-alpha*u))<0){
							gamma = -pow(fabs(((double)1/(2-alpha*u))),((double)1/(n_c+1)));
						}
						else{
							gamma = pow(((double)1/(2-alpha*u)),((double)1/(n_c+1)));
						}
					}
					if(gamma!=gamma || isinf(gamma)){
						i--;
						continue;
						cout<<"FAILURE gamma"<<endl;
						cout<<gamma<<" "<<alpha<<" "<<u<<" "<<n_c<<endl;
						exit(EXIT_FAILURE);
					}
					//cout<<"gamma "<<gamma<<endl;
					double y1 = 0.5*((x1+x2) - gamma*fabs((x2-x1)));
					double y2 = 0.5*((x1+x2) + gamma*fabs((x2-x1)));
					/* cout<<"y1 "<<y1<<endl;
					cout<<"y2 "<<y2<<endl; */
					population_crossover[itr][t][i]=y1;
					population_crossover[itr+1][t][i]=y2;
					if(population_crossover[itr][t][i]!=population_crossover[itr][t][i] || isinf(population_crossover[itr][t][i])){
						population_crossover[itr][t][i] = Q_hmax[i];
						//cout<<x1<<" "<<x2<<" "<<gamma<<endl;
						//cout<<population_crossover[itr+1][t][i]<<endl;
						//cout<<"FAILURE sbx N_h end0"<<endl;
						//exit(EXIT_FAILURE);
					}
					if(population_crossover[itr+1][t][i]!=population_crossover[itr+1][t][i] || isinf(population_crossover[itr+1][t][i])){
						population_crossover[itr+1][t][i] = Q_hmax[i];
						//cout<<x1<<" "<<x2<<" "<<gamma<<endl;
						//cout<<population_crossover[itr+1][t][i]<<endl;
						//cout<<"FAILURE sbx N_h end1"<<endl;
						//exit(EXIT_FAILURE);
					}
				}
			}
			
			for(int i=N_h;i<N_h+N_s;i++){
				double x1;
				double x2;
				double min;
				if(population_crossover[itr][t][i]>=population_crossover[itr+1][t][i]){
					x1 = population_crossover[itr][t][i];
					x2 = population_crossover[itr+1][t][i];
				}
				if(population_crossover[itr][t][i]<population_crossover[itr+1][t][i]){
					x1 = population_crossover[itr][t][i];
					x2 = population_crossover[itr+1][t][i];
				}
				if(population_crossover[itr][t][i]!=population_crossover[itr][t][i] || isinf(population_crossover[itr][t][i])){
					population_crossover[itr][t][i] = P_smin;
					//cout<<population_crossover[itr][t][i]<<endl;
					//cout<<"FAILURE sbx N_h 0"<<endl;
					//exit(EXIT_FAILURE);
				}
				if(population_crossover[itr+1][t][i]!=population_crossover[itr+1][t][i] || isinf(population_crossover[itr+1][t][i])){
					population_crossover[itr+1][t][i] = P_smin;
					//cout<<population_crossover[itr][t][i]<<endl;
					//cout<<"FAILURE sbx N_h 1"<<endl;
					//exit(EXIT_FAILURE);
				}
				if(x1==x2){
					continue;
				}
				else{
					if((x1-P_smin)<(x2-P_smax)){
						min = x1-P_smin;
					}
					else{
						min = x2-P_smax;
					}
					double beta = 1 + ((double)2/(x2-x1))*min;
					/* if(beta!=beta){
						cout<<x2<<" "<<x1<<endl;
						cout<<"FAILURE!"<<endl;
						exit(EXIT_FAILURE);
					} */
					//cout<<-(n_c+1)<<endl;
					double asdf = -(n_c+1);
					//double alpha = 2 - pow(beta,asdf);
					double alpha;
					if(beta<0){
						alpha = 2- -pow(fabs(beta),asdf);
					}
					else{
						alpha = 2- pow(fabs(beta),asdf);
					}
					/* cout<<pow(beta,-(n_c+1))<<endl;
					cout<<"n_c "<<n_c<<endl;
					cout<<"beta "<<beta<<endl;
					cout<<"alpha "<<alpha<<endl; */
					double gamma = 0;
					if(u<=(1/alpha)){
						if((alpha*u)<0){
							gamma = -pow(fabs((alpha*u)),((double)1/(n_c+1)));
						}
						else{
							gamma = pow((alpha*u),((double)1/(n_c+1)));
						}
					}
					else{
						if((1/(2-alpha*u))<0){
							gamma = -pow(fabs((double)(1/(2-alpha*u))),((double)1/(n_c+1)));
						}
						else{
							gamma = pow(((double)1/(2-alpha*u)),((double)1/(n_c+1)));
						}
					}
					//cout<<"gamma "<<gamma<<endl;
					double y1 = 0.5*((x1+x2) - gamma*fabs((x2-x1)));
					double y2 = 0.5*((x1+x2) + gamma*fabs((x2-x1)));
					/* cout<<"y1 "<<y1<<endl;
					cout<<"y2 "<<y2<<endl; */
					population_crossover[itr][t][i]=y1;
					population_crossover[itr+1][t][i]=y2;
					if(population_crossover[itr][t][i]!=population_crossover[itr][t][i] || isinf(population_crossover[itr][t][i])){
						population_crossover[itr][t][i] = P_smin;
						/* cout<<population_crossover[itr][t][i]<<endl;
						cout<<"FAILURE sbx N_s end0"<<endl; */
						//exit(EXIT_FAILURE);
					}
					if(population_crossover[itr+1][t][i]!=population_crossover[itr+1][t][i] || isinf(population_crossover[itr+1][t][i])){
						population_crossover[itr+1][t][i] = P_smin;
						/* cout<<population_crossover[itr+1][t][i]<<endl;
						cout<<"FAILURE sbx N_s end1"<<endl; */
						//exit(EXIT_FAILURE);
					}
				}
			}
		}
		
	}
}

void crossover(){
	population_crossover.clear();
	deque<int> indexes;//list of indexes for crossover, if picked, pick another
	//pick crossover parents
	for(int i=P_e;i<N_p;i++){//index list for N_p-P_e
		indexes.push_back(i);
	}
	for(int i=0;i<P_c;i++){
		int itr1 = rand()%(N_p-P_e) - i;//itr1 ranges from P_e to N_p
		if(itr1<0) itr1=0;
		//cout<<indexes[itr1]<<" ";
		//cout<<itr1<<"\t"<<indexes[itr1]<<endl;
		population_crossover.push_back(population[indexes[itr1]]);
		indexes.erase(indexes.begin()+(itr1));
	}//cout<<endl;
	//cout<<"check!"<<endl;
	/* while(itr<P_c){
		int itr1 = rand()%N_p + N_e;//index of candidate for crossover
			for(int i=0;i<itr;i++){
				if(itr1==indexes[i]) continue;
			}
			indexes[itr] = itr1;//store index
			population_crossover.push_back(population[itr1]);
			itr++;
	} */
	
	//SBX
	sbx();
}

void mutation(){
	population_mutation.clear();
	deque<int> indexes;//list of indexes for crossover, if picked, pick another
	//pick crossover parents
	for(int i=P_e;i<N_p;i++){//index list for N_p-P_e
		indexes.push_back(i);
	}
	for(int i=0;i<P_m;i++){
		int itr1 = rand()%(N_p-P_e) - i;//itr1 ranges from P_e to N_p
		if(itr1<0) itr1=0;
		//cout<<itr1<<"\t"<<indexes[itr1]<<endl;
		population_mutation.push_back(population[indexes[itr1]]);
		indexes.erase(indexes.begin()+(itr1));
	}
	//cout<<"check!"<<endl;
	/* while(itr<P_m){
		int itr1 = rand()%100;
		
		if(itr1<P_e){
			continue;
		}
		else{
			for(int i=0;i<itr1;i++){
				if(itr1==indexes[i]) continue;
			}
			indexes[itr] = itr1;
			population_mutation.push_back(population[itr1]);
			itr++;
		}
		
	} */
	for(int i=0;i<P_m;i++){
		for(int t=0;t<T;t++){
			for(int j=0;j<N_h;j++){
				bool teta = rand()%2;
				//cout<<"teta = "<<teta<<endl;
				if(teta==0){
					double ra;
					do{
						ra = (double)rand()/rand();
						ra -= (int) ra;
					}while(isinf(ra) || ra!=ra);
					
					double temporary;
					//cout<<population_mutation[i][t][j]<<endl;
					if(((double)1-((double)t/T))<0){
						temporary = (double)t/T;
						temporary = 1-temporary;
						temporary = -pow(fabs(temporary), b_mutation);
						temporary = 1-pow(ra,temporary);
						population_mutation[i][t][j] = population_mutation[i][t][j] + (Q_hmax[j]-population_mutation[i][t][j]) * temporary;
						//cout<<"temporary0 "<<temporary<<endl;
					}
					else{
						temporary = (double)t/T;
						temporary = 1-temporary;
						temporary = pow(temporary, b_mutation);
						temporary = 1-pow(ra,temporary);
						population_mutation[i][t][j] = population_mutation[i][t][j] + (Q_hmax[j]-population_mutation[i][t][j]) * temporary;
						//cout<<"temporary1 "<<temporary<<endl;
					}
					
					if(population_mutation[i][t][j]!=population_mutation[i][t][j] || isinf(population_mutation[i][t][j])){
						population_mutation[i][t][j] = Q_hmax[j];
						/* cout<<population_mutation[i][t][j]<<endl;
						cout<<"popow0 "<<pow(ra,pow((1-t/T),b_mutation))<<" "<< pow((1-t/T),b_mutation)<<" ra "<<ra <<endl;
						cout<<"FAILURE MUTATION N_h teta0"<<endl; */
						//exit(EXIT_FAILURE);
					}
				}
				if(teta==1){
					double ra;
					do{
						ra = (double)rand()/rand();
						ra -= (int) ra;
					}while(isinf(ra) || ra!=ra);
				
					double temporary;
					//cout<<population_mutation[i][t][j]<<endl;
					if(((double)1-((double)t/T))<0){
						temporary = (double)t/T;
						temporary = 1-temporary;
						temporary = -pow(fabs(temporary), b_mutation);
						temporary = 1-pow(ra,temporary);
						population_mutation[i][t][j] = population_mutation[i][t][j] + (population_mutation[i][t][j] - Q_hmin[j]) * temporary;
						//cout<<"temporary0 "<<temporary<<endl;
					}
					else{
						temporary = (double)t/T;
						temporary = 1-temporary;
						temporary = pow(temporary, b_mutation);
						temporary = 1-pow(ra,temporary);
						population_mutation[i][t][j] = population_mutation[i][t][j] + (population_mutation[i][t][j] - Q_hmin[j]) * temporary;
						//cout<<"temporary1 "<<temporary<<endl;cout<<"CHECK"<<endl;
					}
					//cout<<"CHECK1"<<endl;
					if(population_mutation[i][t][j]!=population_mutation[i][t][j] || isinf(population_mutation[i][t][j])){
						population_mutation[i][t][j] = Q_hmax[j];
						/* cout<<population_mutation[i][t][j]<<endl;
						cout<<"popow1 "<<pow(ra,pow((1-t/T),b_mutation))<<" "<< pow((1-t/T),b_mutation)<<" ra "<<ra <<endl;
						cout<<"FAILURE MUTATION N_h teta1"<<endl; */
						//exit(EXIT_FAILURE);
					}
					//cout<<"CHECK2"<<endl;
				}
			}
			for(int j=N_h;j<N_h+N_s;j++){
				bool teta = rand()%2;
				if(teta==0){
					double ra;
					do{
						ra = (double)rand()/rand();
						ra -= (int) ra;
					}while(isinf(ra) || ra!=ra);
					
					double temporary;
					//cout<<population_mutation[i][t][j]<<endl;
					if(((double)1-((double)t/T))<0){
						temporary = (double)t/T;
						temporary = 1-temporary;
						temporary = -pow(fabs(temporary), b_mutation);
						temporary = 1-pow(ra,temporary);
						population_mutation[i][t][j] = population_mutation[i][t][j] + (P_smax-population_mutation[i][t][j]) * temporary;
						//cout<<"temporary0 "<<temporary<<endl;
					}
					else{
						temporary = (double)t/T;
						temporary = 1-temporary;
						temporary = pow(temporary, b_mutation);
						temporary = 1-pow(ra,temporary);
						population_mutation[i][t][j] = population_mutation[i][t][j] + (P_smax-population_mutation[i][t][j]) * temporary;
						//cout<<"temporary1 "<<temporary<<endl;
					}
					if(population_mutation[i][t][j]!=population_mutation[i][t][j] || isinf(population_mutation[i][t][j])){
						population_mutation[i][t][j] = P_smin;
						/* cout<<population_mutation[i][t][j]<<endl;
						cout<<"popow0 "<<pow(ra,pow((1-t/T),b_mutation))<<" "<< pow((1-t/T),b_mutation)<<" ra "<<ra <<endl;
						cout<<"FAILURE MUTATION N_s teta0"<<endl; */
						//exit(EXIT_FAILURE);
					}
				}//cout<<"CHECK3"<<endl;
				if(teta==1){
					double ra;
					do{
						ra = (double)rand()/rand();
						ra -= (int) ra;
					}while(isinf(ra) || ra!=ra);
					
					//cout<<population_mutation[i][t][j]<<endl;
					double temporary;
					if(((double)1-((double)t/T))<0){
						temporary = (double)t/T;
						temporary = 1-temporary;
						temporary = -pow(fabs(temporary), b_mutation);
						temporary = 1-pow(ra,temporary);
						population_mutation[i][t][j] = population_mutation[i][t][j] + (population_mutation[i][t][j]-P_smin) * temporary;
						//cout<<"temporary0 "<<temporary<<endl;
					}
					else{
						temporary = (double)t/T;
						temporary = 1-temporary;
						temporary = pow(temporary, b_mutation);
						temporary = 1-pow(ra,temporary);
						population_mutation[i][t][j] = population_mutation[i][t][j] + (population_mutation[i][t][j]-P_smin) * temporary;
						//cout<<"temporary1 "<<temporary<<endl;
					}
					if(population_mutation[i][t][j]!=population_mutation[i][t][j] || isinf(population_mutation[i][t][j])){
						population_mutation[i][t][j] = P_smin;
						/* cout<<population_mutation[i][t][j]<<endl;
						cout<<"popow0 "<<pow(ra,pow((1-t/T),b_mutation))<<" "<< pow((1-t/T),b_mutation)<<" ra "<<ra <<endl;
						cout<<"FAILURE MUTATION N_s teta1"<<endl; */
						//exit(EXIT_FAILURE);
					}
				}//cout<<"CHECK4"<<endl;
			}
		}
	}
}

void update_population(){
	for(int i=0;i<P_e;i++){//add elite population to P
		population[i].swap(population_elite[i]);
	}//cout<<population_crossover.size()<<" ";
	for(int i=0;i<P_c;i++){//add crossover population to P
		population[P_e+i].swap(population_crossover[i]);
	}//cout<<"huehue"<<endl;
	for(int i=0;i<P_m;i++){
		population[i+P_e+P_c].swap(population_mutation[i]);
	}
	
}

double evaluation(deque<deque<double> > temp1){
	double temp=0;
	for(int t=0;t<T;t++){
		for(int i=0;i<N_h+N_s;i++){
			temp = temp + a_si + b_si*temp1[t][N_h+i] + c_si*pow(temp1[t][N_h+i],2);// + fabs(d_si*sin(e_si*(P_smin - temp1[t][N_h+i])));
		}
	}
	return temp;
}

void evaluation_pso(){//calculate fitness value for population
	//cout<<"evaluation_pso"<<endl;
	population_pso_fitness.clear();
	for(int itr=0;itr<POP_SIZE;itr++){
		double temp=0;
		for(int t=0;t<T;t++){
			for(int i=0;i<N_s;i++){
				//cout<<itr<<" "<<t<<" "<<i<<endl;
				temp = temp + a_si + b_si*population_pso[itr][t][N_h+i] + c_si*pow(population_pso[itr][t][N_h+i],2);// + fabs(d_si*sin(e_si*(P_smin - population[itr][t][N_h+i])));
				//cout<<temp<<" ";
			}
			//cout<<endl;
			
		}

		population_pso_fitness.push_back(temp);//default place for the fitness value
	}
	
	//bubble sort, no sorting
/* 	for(int itr1=0;itr1<N_p;itr1++){
		for(int itr2=0;itr2<N_p;itr2++){
			if(population_pso_fitness[itr1]<population_pso_fitness[itr2]){
				population_pso[itr1].swap(population_pso[itr2]);
				double huehue = population_pso_fitness[itr1];
				population_pso_fitness[itr1] = population_pso_fitness[itr2];
				population_pso_fitness[itr2] = huehue;
				//population_fitness[itr1].swap(population_fitness[itr2]);
			}
		}
	} */
	//cout<<"evaluation pso end"<<endl;
}

void initialize_pso(){//initialize all velocity to 0
	//cout<<"initialize_pso"<<endl;
//initialize initial position as pbest(list of best position)
	for(int i=0;i<N_p;i++){
		deque<deque<double> > temp1;
		for(int t=0;t<T;t++){
			deque<double> temp2;
			
			for(int j=0;j<N_h;j++){
				bool ran = rand()%2;
				double temp3 = abs(static_cast<double> (rand() / static_cast<double> (RAND_MAX/(Q_hmin[j] - Q_hmax[j]) )));
				//double temp3 = abs(rand()%(Q_hmin[j] - Q_hmax[j]));
				if(ran==true){
					temp3 = -temp3;
				}
				if(isinf(temp3) || temp3!=temp3){
					j--;
					continue;
					cout<<"FAILURE INIT TEMP3"<<endl;
					cout<<temp3<<endl;
					exit(EXIT_FAILURE);
				}	
				temp2.push_back(temp3);
			}
			
			for(int j=N_h;j<N_h+N_s;j++){
				bool ran = rand()%2;
				double temp3 = abs(static_cast<double> (rand() / static_cast<double> (RAND_MAX/(Q_hmin[j] - Q_hmax[j]) )));
				if(ran==true){
					temp3 = -temp3;
				}
				if(isinf(temp3) || temp3!=temp3){
					j--;
					continue;
					cout<<"FAILURE INIT TEMP3"<<endl;
					cout<<temp3<<endl;
					exit(EXIT_FAILURE);
				}	
				temp2.push_back(temp3);
			}
			
			temp1.push_back(temp2);
			temp2.clear();
		}
			velocity.push_back(temp1);
			temp1.clear();
			//pbest[i] = population[i];
	}
	
	for(int i=0;i<N_p;i++){//copy population to pbest
		deque<deque<double> > temp1;
		for(int t=0;t<T;t++){
			deque<double> temp2;
			for(int j=0;j<N_h+N_s;j++){
				temp2.push_back(population[i][t][j]);
			}
			temp1.push_back(temp2);
			temp2.clear();
		}
		pbest.push_back(temp1);
		temp1.clear();
		pbest_fitness.push_back(population_fitness[i]);
	}
		//initialize population_pso to zero
	for(int i=0;i<N_p;i++){
		deque<deque<double> > temp1;
		for(int t=0;t<T;t++){
			deque<double> temp2;
			for(int j=0;j<N_h+N_s;j++){
				temp2.push_back(0);
			}
			temp1.push_back(temp2);
			temp2.clear();
		}
		population_pso.push_back(temp1);
		temp1.clear();
	}
		//get first as gbest
		gbest = population[0];
		gbest_fitness = population_fitness[0];
		//cout<<"gbest size"<<gbest.size()<<endl;
		//cout<<"initialize pso end"<<endl;
}

double update_velocity(int i, int t, int j, int itr){//returns updated velocity of i t j
		//cout<<"update_velocity"<<endl;
		double c1, c2;
		double c1i = 2.5;
		double c1f = 0.5;
		double c2i = 0.5;
		double c2f = 2.5;
		double rand1, rand2, rand3;
		double newv;
		do{
			rand1 = (double)rand()/rand();
			rand1 -= (int)rand1;
		}while(isinf(rand1) || rand1!=rand1);
		do{
			rand2 = (double)rand()/rand();
			rand2 -= (int)rand2;
		}while(isinf(rand2) || rand2!=rand2);
		do{
			rand3 = (double)rand()/rand();
			rand3 -= (int)rand3;
		}while(isinf(rand3) || rand3!=rand3);
		
		if(velocity[i][t][j]==0 && rand3<0.5){
			double rand4;
			do{
			rand4 = (double)rand()/rand();
			rand4 -= (int)rand4;
			}while(isinf(rand4) || rand4!=rand4);
			if(j<4){//vmax is for hydro
				newv = rand4*Q_hmax[j];
			}
			else{//vmax is for thermal
				newv = rand4*P_smax;
			}
			return newv;
		}
		else{
			double rand5;
			do{
			rand5 = (double)rand()/rand();
			rand5 -= (int)rand5;
			}while(isinf(rand5) || rand5!=rand5);
			if(j<4){//vmax is for hydro
				newv = -rand5*Q_hmax[j];
			}
			else{//vmax is for thermal
				newv = -rand5*P_smax;
			}
			return newv;
		}
		c1 = (c1f-c1i)*(double)(itr/gmax_pso) + c1i;
		c2 = (c2f-c2i)*(double)(itr/gmax_pso) + c2i;
		newv = c1*rand1*(pbest[i][t][j] - population[i][t][j]) + c2*rand2*(gbest[t][j] - population[i][t][j]);
		
		
/* 		if(isinf(C) || C!=C){
			cout<<"FAILURE C"<<endl;
			cout<<C<<endl;
			exit(EXIT_FAILURE);
		}
		if(isinf(w) || w!=w){
			cout<<"FAILURE w"<<endl;
			cout<<w<<endl;
			exit(EXIT_FAILURE);
		} */
		//cout<<"res: "<<pbest[i][t][j]<<endl;
		//cout<<"check"<<endl;
/* 		if(isinf(newv) || newv!=newv){
			cout<<"FAILURE OUTPUT"<<endl;
			cout<<"velocity: "<<velocity[i][t][j]<<endl;
			cout<<"pbest: "<<pbest[i][t][j]<<endl;
			cout<<"population: "<<population[i][t][j]<<endl;
			cout<<"gbest: "<<gbest[t][j]<<endl;
			cout<<newv<<endl;
			exit(EXIT_FAILURE);
		} */
		//cout<<"update velocity end"<<endl;
		return newv;		
}

void pso(int itr){//run pso for all population
	//cout<<"pso"<<endl;
//change velocity deque	
	for(int i=0;i<N_p;i++){
		for(int t=0;t<T;t++){
			for(int j=0;j<N_h+N_s;j++){
				velocity[i][t][j] = update_velocity(i,t,j,itr);
				population_pso[i][t][j] = population[i][t][j] + velocity[i][t][j];
	/* 			if(isinf(population_pso[i][t][j]) || population_pso[i][t][j]!=population_pso[i][t][j]){
					cout<<"FAIL PSO main pop"<<endl;
					cout<<population_pso[i][t][j]<<endl;
					exit(EXIT_FAILURE);
				}
				if(isinf(velocity[i][t][j]) || velocity[i][t][j]!=velocity[i][t][j]){
					cout<<"FAIL PSO main vel"<<endl;
					cout<<velocity[i][t][j]<<endl;
					exit(EXIT_FAILURE);
				} */
				//cout<<"check1"<<endl;
			}
			//cout<<"check2"<<endl;
		}
	}
	//evaluation
	evaluation_pso();
	
	//update pbest
	for(int i=0;i<N_p;i++){
		if(population_fitness[i]<pbest_fitness[i]){
			for(int t=0;t<T;t++){
				for(int j=0;j<N_h+N_s;j++){
					pbest[i][t][j] = population[i][t][j];
				}
			}
			pbest_fitness[i] = population_fitness[i];
		}
		else{
			continue;
		}	
	}
	
	//update population if necessary from population_pso
	for(int i=0;i<N_p;i++){
		if(population_pso_fitness[i] < population_fitness[i]){
			population[i].swap(population_pso[i]);
			population_fitness[i] = population_pso_fitness[i];
		}
	}
	
	//sort all individuals, no sorting
/* 	for(int itr1=0;itr1<POP_SIZE;itr1++){
		for(int itr2=0;itr2<POP_SIZE;itr2++){
			if(population_fitness[itr1]<population_fitness[itr2]){
				population[itr1].swap(population[itr2]);
				double huehue = population_fitness[itr1];
				population_fitness[itr1] = population_fitness[itr2];
				population_fitness[itr2] = huehue;
				//population_fitness[itr1].swap(population_fitness[itr2]);
			}
		}
	} */
	//update gbest
	double bestIndex;
	for(int i=0;i<N_p;i++){
		double temp = 999999999999999;
		if(population_fitness[i]<temp){
			temp = population_fitness[i];
			bestIndex = i;
		}
	}
	gbest_fitness = population_fitness[bestIndex];
	for(int t=0;t<T;t++){
		for(int j=0;j<N_h+N_s;j++){
			gbest[t][j] = population[bestIndex][t][j];
		}
	}
	//cout<<"pso end"<<endl;
}

void print_pop(){
	for(int i =0;i<N_p;i++){
		for(int t=0;t<T;t++){
			for(int j=0;j<N_h+N_s;j++){
				cout<<population[i][t][j]<<" ";
			}
			cout<<endl;
		}
		cout<<endl;
	}	
}

void constraints_handling(){
//cout<<"constraints handling"<<endl;
	for(int i=0;i<N_p;i++){
		for(int j=0;j<N_h;j++){//constraints handling for hydro plants
			int iteration;
			iteration = 0;
			check_discharge(i,j);//check Q_h(j,t)
			double V_hc = volume[i][T-1][j] - V_hend[j];
			//cout<<"V_hc = "<< V_hc<<endl;
			while(fabs(V_hc)>e_vcourse && iteration<=iteration_course){
				double V_haverage = V_hc/T;
				//cout<<"V_haverage = "<<V_haverage<<endl;
				for(int t=0;t<T;t++){
					population[i][t][j] = population[i][t][j] + V_haverage;
					//cout<<"population[][][]: "<<population[i][t][j]<<endl;
					if(population[i][t][j]!=population[i][t][j] || isinf(population[i][t][j])){
						population[i][t][j] = Q_hmax[j];
						cout<<"FAILURE N_h course"<<endl;
						exit(EXIT_FAILURE);
					}
					check_discharge1(i,j,t);
				}
				//cout<<"start update volume i 1"<<endl;
				update_volume(i);
				//cout<<"end update volume i 1"<<endl; 
				V_hc = volume[i][T-1][j] - V_hend[j];
				//cout<<"V_hc = "<< V_hc<<endl;
				iteration++;
			}
			//cout<<"end while course"<<endl;
			//cout<<"V_hc: "<<V_hc<<endl;
			//cout<<"checkC "<< i <<endl;
			iteration=0;
			deque<double> Q_hc;
			while(fabs(V_hc)>e_vfine && iteration<iteration_fine){
				//cout<<"iteration: "<<iteration<<endl;
				Q_hc.clear();
				Q_hc = Q_hcF(i,j);//get the change in water discharge population from i and j
				int itr = max_adjustableQ(Q_hc);
				//cout<<"itr: "<<itr<<endl;
				//cout<<"population[i][itr][j] "<<population[i][itr][j]<<endl;
				population[i][itr][j] = population[i][itr][j] + V_hc;
				//cout<<"V_hc: "<<V_hc<<endl;
				//cout<<"pop: "<<population[i][itr][j]<<endl;
				if(population[i][itr][j]!=population[i][itr][j] || isinf(population[i][itr][j])){
					population[i][itr][j] = Q_hmax[j];
					cout<<itr<<endl;
					cout<<population[i][itr][j]<< " "<<V_hc<<endl;
					cout<<"FAILURE N_h fine"<<endl;
					exit(EXIT_FAILURE);
				}
				//cout<<"start update volume i 2"<<endl;
				update_volume(i);
				//cout<<"end update volume i 2"<<endl;
				check_discharge1(i,j,itr);
				//cout<<"end check discharge1 "<<endl;
				V_hc = volume[i][T-1][j] - V_hend[j];
				//cout<<"end V_hc"<<endl;
				iteration++;
			}
			//cout<<"end while fine"<<endl;
		}
		//cout<<"end constraints hydro "<<endl;
		//cout<<"CONSTRAINTS FUNCTION"<<endl;
		//cout<<"start update volume i 3"<<endl;
		update_volume(i);
		//cout<<"end update volume i 3"<<endl;
		//cout<<"update volume fine end"<<endl;
		for(int t=0;t<T;t++){//thermal plants constraints handling
			int iteration = 0;
			double P_sc = P_D[t] - P_ssum(i,t) - P_hsum(i,t);//no transmission loss
			population[i][t][4] += P_sc;
			/* while(fabs(P_sc)>e_pcourse && iteration<=iteration_course){
				double P_taverage = P_sc/N_s;
				
				for(int ns = N_h; ns<N_h+N_s;ns++){
					population[i][t][ns] = population[i][t][ns] + P_taverage;
					if(population[i][t][ns]!=population[i][t][ns] || isinf(population[i][t][ns])){
						population[i][t][ns] = P_smin;
						//cout<<"FAILURE N_s course"<<endl;
						//exit(EXIT_FAILURE);
					}
					check_thermal1(i,t,ns);
				}
				P_sc = P_D[t] - P_ssum(i,t) - P_hsum(i,t);//no transmission loss
				iteration++;
			}
			iteration =0;
			//deque<double> P_sc1;
			while(fabs(P_sc)>e_pfine && iteration<=iteration_fine){
				//P_sc1.clear();
				P_sc = P_D[t] - P_ssum(i,t) + P_hsum(i,t);
				//P_sc1 = P_scF(i,t);//i = pop_index, t = time
				int itr = 4;//there is only one thermal generation that will make up for the power deficiency
				//cout<<"ITR: "<<itr<<endl;
				population[i][t][itr] = population[i][t][itr] + P_sc;
				if(population[i][t][itr]!=population[i][t][itr] || isinf(population[i][t][itr])){
					population[i][t][itr] = P_smin;
					//cout<<"FAILURE N_s fine"<<endl;
					//exit(EXIT_FAILURE);
				}
				P_sc = P_D[t] - P_ssum(i,t) - P_hsum(i,t);
				iteration++;
			} */
			//P_sc1.clear();
		}
		//cout<<"end constraints thermal"<<endl;
	}
	//cout<<"CONSTRAINTS END"<<endl;
	for(int i=0;i<POP_SIZE;i++){
		for(int t=0;t<T;t++){
			for(int j=0;j<N_h;j++){
				if(isinf(population[i][t][j]) || (population[i][t][j]!=population[i][t][j])){
					population[i][t][j] = Q_hmax[j];
					cout<<"FAILURE C3"<<endl;
					exit(EXIT_FAILURE);
				}
			}
			for(int j=N_h;j<N_h+N_s;j++){
				if(isinf(population[i][t][j]) || (population[i][t][j]!=population[i][t][j])){
					population[i][t][j] = P_smin;
					cout<<"FAILURE C3"<<endl;
					exit(EXIT_FAILURE);
				}
			}
		}
	}
	//cout<<"constraints handling end"<<endl;
}


int main(){
best_fitness = 999999999999999;
inflowF();
initUpstream();
aveTime = 0;
for(int t=0;t<T;t++){//initialization of super global variable
	deque<double> temps;//chromosome
	deque<double> temps1;//Ph Qh
	deque<double> temps2;//Ps
	if(t==0){	
		for(int j=0;j<N_h+N_s;j++){
			temps.push_back(0);
		}
		for(int j=0;j<N_h;j++){
			temps1.push_back(0);
		}
		for(int j=N_h;j<N_h+N_s;j++){
			temps2.push_back(0);
		}
	}
	aveBestChromosome.push_back(temps);
	aveWorstChromosome.push_back(temps);
	aveBestPh.push_back(temps1);
	aveBestPs.push_back(temps2);
	aveBestQh.push_back(temps1);
}



for(int mainItr=0;mainItr<iterators;mainItr++){
	population.clear();
	population_fitness.clear();
	volume.clear();
	population_elite.clear();
	population_crossover.clear();
	population_mutation.clear();
	velocity.clear();
	population_pso.clear();
	population_pso_fitness.clear();
	pbest.clear();
	gbest.clear();
	gbest_fitness = 999999999999999;
	best_fitness = 9999999999999999;
	srand(time(NULL));

	//INITIALIZATION 
	for(int i=0;i<POP_SIZE;i++){//temp(n), n - dimension
		deque<deque<double> > temp2;
		for(int t=0;t<T;t++){
			deque<double> temp1;
			
			for(int nh=0;nh<N_h;nh++){
				double ran;
				double temp;
				do{
					ran = (double)rand()/rand();
					ran -= (int) ran;
				}while(isinf(ran) || ran!=ran);
				temp = Q_hmin[nh] + ran*(Q_hmax[nh]-Q_hmin[nh]);
				//cout<<temp<<" ";
				temp1.push_back(temp);
			}
			
			for(int ns=N_h;ns<N_h+N_s;ns++){
				double ran;
				double temp;
				do{
					ran = (double)rand()/rand();
					ran -= (int) ran;
				}while(isinf(ran) || ran!=ran);
				temp = P_smin + ran*(P_smax-P_smin);
				//cout<<temp<<" ";
				temp1.push_back(temp);
			} temp2.push_back(temp1);
			temp1.clear();
		}
		population.push_back(temp2);
		temp2.clear();
	}


	for(int i=0;i<POP_SIZE;i++){//water dynamic volume generation
		//water dynamic balance generation
		deque<deque<double> > V_h;//V_h(t,j)
		V_h.clear();
		for(int t=0;t<T;t++){
			deque<double> temp1;
			for(int j=0;j<N_h;j++){
				double temp=0;//V_h(j)
				if((t-1)<0){
					temp = V_hbegin[j];// + inflow[j][t] - population[i][t][j] + sum_upstream(i,t,j);//spillage not calculated
					//cout<<i<<" "<<t<<" "<<j<<endl;
					//cout<<temp<<"\tV_hbegin: "<<V_hbegin[j]<<endl;//" inflow: "<<inflow[j][t]<<" discharge: "<<population[i][t][j]<<" upstream: "<<sum_upstream(i,t,j)<<endl;
					temp1.push_back(temp);
				}
				else{
					temp = V_h[t-1][j] + inflow[j][t] - population[i][t][j] + sum_upstream(i,t,j);//spillage not calculated
					//cout<<i<<" "<<t<<" "<<j<<endl;
					//cout<<temp<<"\tV[t-1]: "<<V_h[t-1][j]<<" inflow: "<<inflow[j][t]<<" discharge: "<<population[i][t][j]<<" upstream: "<<sum_upstream(i,t,j)<<endl;
					//if(temp<0) cout<<"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"<<endl;
					temp1.push_back(temp);
				}
			}
			V_h.push_back(temp1);
			temp1.clear();
		}
		volume.push_back(V_h);//always update volume reservoir after change in values
		V_h.clear();
	}

	//CONSTRAINTS HANDLING
	for(int i=0;i<POP_SIZE;i++){
		for(int j=0;j<N_h;j++){//constraints handling for hydro plants
			int iteration;
			iteration = 0;
			check_discharge(i,j);//check Q_h(j,t)
			double V_hc = volume[i][T-1][j] - V_hend[j];
			//cout<<"V_hc = "<< V_hc<<endl;
			while(fabs(V_hc)>e_vcourse && iteration<iteration_course){
				double V_haverage = V_hc/T;
				//cout<<"V_haverage = "<<V_haverage<<endl;
				for(int t=0;t<T;t++){
					population[i][t][j] = population[i][t][j] + V_haverage;
					/* if(population[i][t][j]!=population[i][t][j]){
							cout<<"main FAILURE N_h course"<<endl;
							exit(EXIT_FAILURE);
					} */
					check_discharge1(i,j,t);
				}
				update_volume(i);
				V_hc = volume[i][T-1][j] - V_hend[j];
				iteration++;
			}
			iteration=0;
			deque<double> Q_hc;
			while(fabs(V_hc)>e_vfine && iteration<iteration_fine){
				Q_hc.clear();
				Q_hc = Q_hcF(i,j);//get the population from i and j, change in Q
				int itr = max_adjustableQ(Q_hc);
				population[i][itr][j] = population[i][itr][j] + V_hc;
				/* if(population[i][itr][j]!=population[i][itr][j]){
						cout<<itr<<endl;
						cout<<population[i][itr][j]<< " "<<V_hc<<endl;
						cout<<"main FAILURE N_h fine"<<endl;
						exit(EXIT_FAILURE);
				} */
				update_volume(i);
				check_discharge1(i,j,itr);
				V_hc = volume[i][T-1][j] - V_hend[j];
				iteration++;
			}
		}
		update_volume(i);
		
		for(int t=0;t<T;t++){//thermal plants constraints handling
			int iteration = 0;
			double P_sc = P_D[t] - P_ssum(i,t) - P_hsum(i,t);//no transmission loss
			population[i][t][4] +=P_sc;
			/* while(fabs(P_sc)>e_pcourse && iteration<=iteration_course){
				double P_taverage = P_sc/N_s;
				
				for(int ns = N_h; ns<N_s;ns++){
					population[i][t][ns] = population[i][t][ns] + P_taverage;
					check_thermal1(i,t,ns);
				}
				P_sc = P_D[t] - P_ssum(i,t) - P_hsum(i,t);//no transmission loss
				iteration++;
			}
			iteration =0;
			//deque<double> P_sc1;//change in thermal plants power 
			while(fabs(P_sc)>e_pfine && iteration<=iteration_fine){
				//P_sc1.clear();
				P_sc = P_D[t] - P_ssum(i,t) - P_hsum(i,t);//violation in the equality constraints of power generated to power demand
				//P_sc1 = P_scF(i,t);//i = pop_index, t = time
				int itr = 4;//there is only one thermal plant
				population[i][t][itr] = population[i][t][itr] + P_sc;
				P_sc = P_D[t] - P_ssum(i,t) - P_hsum(i,t);
				iteration++;
			} */
			//P_sc1.clear();
		}
	}
	clock_t t;
	t = clock();
	//RCGA
	for(int g = 0;g<gmax_rcga;g++){
		//EVALUATION
		evaluation();//generate fitness value and sort them
		
		elitist();//copy all the elite population, get best_fitness

		//CROSSOVER
		crossover();//pure crossover w/ sbx. 

		//MUTATION
		mutation();//pure mutation

		update_population();//add the crossover, mutated population to the general population
		
		constraints_handling();
		
		population_fitness.clear();
		population_elite.clear();
		population_crossover.clear();
		population_mutation.clear();
	}
	constraints_handling();

	cout<<"MAINITR: "<<mainItr<<endl;
	//PSO
	initialize_pso();
	for(int g = 0;g<gmax_pso;g++){
	//cout<<"gpso: "<<g<<endl;
		pso(g);
		//constraints handling for pso
		constraints_handling();
		//print_pop();
	}
	t=clock() - t;

	evaluation();
	cout<<"----------------------------------------------------------------"<<endl;
	aveBestFitness.push_back(population_fitness[0]);
	aveWorstFitness.push_back(population_fitness[99]);
	
	cout<<"Best Fitness: "<<population_fitness[0]<<endl;
	cout<<"Worst Fitness: "<<population_fitness[99]<<endl;
	
	for(int t=0;t<T;t++){
		for(int j=0;j<N_h+N_s;j++){
			aveBestChromosome[t][j]+=population[0][t][j];
		}
	}
	for(int t=0;t<T;t++){
		for(int j=0;j<N_h+N_s;j++){
			aveWorstChromosome[t][j]+=population[99][t][j];
		}
	}
	cout<<"----------------------------------------------------------------"<<endl;
	cout<<"P_h: "<<endl;
	for(int t=0;t<T;t++){//add average best P_h
		for(int a=0;a<N_h;a++){
			double temp = 0;
			temp = c1[a]*pow(volume[0][t][a],2) + c2[a]*pow(population[0][t][a],2) + c3[a]*volume[0][t][a]*population[0][t][a] + c4[a]*volume[0][t][a] + c5[a]*population[0][t][a] + c6[a];
			if(temp<0){
				temp=0;
			}cout<<temp<<" ";
			aveBestPh[t][a] += temp;
		}
		cout<<endl;
	}	
	cout<<"----------------------------------------------------------------"<<endl;	
	for(int t=0;t<T;t++){
		for(int j=N_h;j<N_h+N_s;j++){
			aveBestPs[t][j] += population[0][t][j];
		}
	}
	
	for(int t=0;t<T;t++){
		for(int j=0;j<N_h;j++){
			aveBestQh[t][j] += population[0][t][j];
		}
	}
	
	cout<<"Best: "<<endl;
	for(int i=0;i<1;i++){
		for(int t=0;t<T;t++){
			for(int j=0;j<N_h+N_s;j++){
				cout<<population[i][t][j]<<" ";
			}
			cout<<endl;
		}
		cout<<endl;
	}
	cout<<endl;
	
	
	cout<<"Volume: "<<endl;
	for(int i=0;i<2;i++){
		for(int t=0;t<T;t++){
			for(int j=0;j<N_h;j++){
				cout<<volume[i][t][j]<<" ";
			}
			cout<<endl;
		}
		cout<<endl;
	}
	cout<<endl;
	cout<<"----------------------------------------------------------------"<<endl; 
	aveTime += (float)t/CLOCKS_PER_SEC;
	cout<<(float)t/CLOCKS_PER_SEC<<" secs"<<endl;
}//end of main iteration (30 runs)


for(int t=0;t<T;t++){
	for(int j=0;j<N_h+N_s;j++){
		aveBestChromosome[t][j] = (double)aveBestChromosome[t][j]/iterators;
		aveWorstChromosome[t][j] = (double)aveWorstChromosome[t][j]/iterators;
	}
}

for(int t=0;t<T;t++){
	for(int j=0;j<N_h;j++){
		aveBestPh[t][j] = (double)aveBestPh[t][j]/iterators;
		aveBestQh[t][j] = (double)aveBestQh[t][j]/iterators;
	}
}

for(int t=0;t<T;t++){
	for(int j=0;j<N_s;j++){
		aveBestPs[t][j] = (double)aveBestPs[t][j]/iterators;
	}
}

double aveBestF=0;
double aveWorstF=0;
for(int i=0;i<30;i++){
	aveBestF += aveBestFitness[i];
	aveWorstF += aveWorstFitness[i];
}
aveBestF = (double)aveBestF/iterators;
aveWorstF = (double)aveWorstF/iterators;

cout<<"Average Best: "<<aveBestF<<endl;
cout<<"Average Worst: "<<aveWorstF<<endl;
cout<<endl;
cout<<"============================================================"<<endl;
cout<<"Average Best Chromosome: "<<endl;
for(int t=0;t<T;t++){
	for(int j=0;j<N_h+N_s;j++){
		cout<<aveBestChromosome[t][j]<<" ";
	}
	cout<<endl;
}
cout<<"============================================================"<<endl;
cout<<"Average Worst Chromosome: "<<endl;
for(int t=0;t<T;t++){
	for(int j=0;j<N_h+N_s;j++){
		cout<<aveWorstChromosome[t][j]<<" ";
	}
	cout<<endl;
}
cout<<"============================================================"<<endl;
cout<<"Average Best P_h: "<<endl;
for(int t=0;t<T;t++){
	for(int j=0;j<N_h;j++){
		cout<<aveBestPh[t][j]<<" ";
	}
	cout<<endl;
}
cout<<"============================================================"<<endl;
cout<<"Average Best P_s"<<endl;
for(int t=0;t<T;t++){
	for(int j=0;j<N_s;j++){
		cout<<aveBestPs[t][j]<<" ";
	}
	cout<<endl;
}
cout<<"============================================================"<<endl;
cout<<"Total Generation:"<<endl;
for(int t=0;t<T;t++){
	double temp = 0;
	for(int j=0;j<N_h;j++){
		temp += aveBestPh[t][j];
	}
	for(int j=N_h;j<N_h+N_s;j++){
		temp += aveBestChromosome[t][j];
	}
	cout<<temp<<endl;
}
cout<<"============================================================"<<endl;
cout<<"Average Best Q_h"<<endl;
for(int t=0;t<T;t++){
	for(int j=0;j<N_h;j++){
		cout<<aveBestQh[t][j]<<" ";
	}
	cout<<endl;
}
cout<<"============================================================"<<endl;





return 0;
}










