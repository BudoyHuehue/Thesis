#include <iostream>
#include <fstream>
#include <string>
#include <deque>
#include <cstdlib>
#include <cstdio>
#include <cctype>
#include <ctime>
#include <iomanip>
using namespace std;

double delay(int a, int b,deque< deque< int > > fmax, deque< deque < int > > fcurrent){//get the fmax and fcurrent from deque fmax & fcurrent
	double m,c;
	
	m = fmax[a][b];
	c = fcurrent[a][b];
	if((m-c)==0) return -100;
	else return c/(m-c);

}

int getfcurrent(int index1, int index2){

}
int getfmax(int index1, int index2){
	
}
/*double delay_sum(deque< int > chromosomes, deque<deque<int> > fcurrent, deque<deque<int> > fmax,deque<int> &delay){
	double tempint;//
	for(int i=0;i<chromosomes.size()-1;i++){
		deque<int> temp = chromosomes[i];
		for(int j=0;j<temp.size();j++){
			fcurrent[chromosomes[i][j]][chromosomes[i][j+1]]
			fmax[chromosomes[i][j]][chromosomes[i][j+1]]
		}
	cout<<"check"<<endl;
		
		cout<<delay(chromosomes[i],chromosomes[i+1])<<" ";
	}
	cout<<endl;
	//cout<<std::fixed<<std::setprecision(1000)<<tempint;
	return tempint;
}*/


int main (){
srand(time(NULL));
deque< deque <int> > chromosomes;//population
deque< deque <int> > kids;
deque<int> delay_array;
int itr;

ifstream myfile ("fmax.txt");
deque< deque<int> > fmax;
if(myfile.is_open()){
	string line;
	while (getline(myfile,line)){//read line and store
		deque<int> temp;
		for(int i=0;i<line.length();i++){
			string str_temp;
			while(isdigit(line[i])){
				str_temp = str_temp+line[i];
				//temp.push_back((int)line[i]);
				i++;
			}
			temp.push_back(atoi(str_temp.c_str()));
			//cout<<atoi(str_temp.c_str())<<endl;
		}
		fmax.push_back(temp);
	}
	myfile.close();
}
else cout<<"Can't open file"<<endl;

ifstream myfile1 ("fcurrent2.txt");
deque< deque<int> > fcurrent;
if(myfile1.is_open()){
	string line;
	while (getline(myfile1,line)){//read line and store
		deque<int> temp;
		for(int i=0;i<line.length();i++){
			string str_temp;
			while(isdigit(line[i])){
				str_temp = str_temp+line[i];
				//temp.push_back((int)line[i]);
				i++;
			}
			temp.push_back(atoi(str_temp.c_str()));
			//cout<<atoi(str_temp.c_str())<<endl;
		}
		fcurrent.push_back(temp);
	}
	myfile1.close();
}
else cout<<"Can't open file"<<endl;

deque<int> max = fcurrent[0];


int start_point, end_point;
cout<<"Enter start point: ";
cin>>start_point;
cout<<endl;
cout<<"Enter end point: ";
cin>>end_point;
cout<<endl;

//INITIALIZE POPULATION-------------------------------------------------------

for(int i=0;i<40;i++){//initialize population
//cout<<i<<endl;
	bool flag = 0;
	deque < int > chromosomes_temp;
	
	chromosomes_temp.push_back(start_point);
	//cout<<start_point<<" ";
	//checking of blocked paths not done yet
	while(1){//populate the population until end_point shows
		int temp = rand()%fcurrent.size();

		if(end_point!=temp){
			for(int k=0;k<chromosomes_temp.size();k++){//check redundancy of chromosome
				if((temp == chromosomes_temp[k])||(temp==start_point)||(temp==end_point)&&(delay(chromosomes_temp[k-1],chromosomes_temp[k],fmax,fcurrent)<0)){//break if chromosome already exists, raise flag to ignore temp
					//cout<<temp<<"break--redundant";
					k--;
					flag = 1;
					break;
				}
			}
		}
		
		if(end_point==temp){
			//cout<<temp<<" ";
			chromosomes_temp.push_back(temp);
			break;
		}
		
		if(flag==1){
			flag=0;
			continue;
		}
		
		if(flag==0){//end the randomization of chromosome
			//cout<<temp<<" ";
			chromosomes_temp.push_back(temp);
		}
	}
	//cout<<endl;
	chromosomes.push_back(chromosomes_temp);
	chromosomes_temp.clear();
	//cout<<endl;
}

//EVALUATION
/*for(int i=0;i<40;i++){
	cout<<"\n NEW: ["<<i<<"]"<<endl;
	delay_sum(chromosomes[i],fcurrent,fmax,delay);
	cout<<endl;
}*/
/* for(int i=0;i<40;i++){
	deque < int > temp = chromosomes[i];
	for(int j=0;j<temp.size();j++){
		cout<<temp[j];
	}
cout<<endl;
} */
cout<<endl;
cout<<endl;
int y=0;
double maindelaybest=99999;
for(itr=0;itr<1000000;itr++){

maindelaybest=99999;
if(y==999) break;
double delayvaluebest=9999;


for(int i=0;i<40;i++){//calculate for best delay
	deque < int > temp = chromosomes[i];
	
	double delayvalue=0;
	for(int j=0;j<temp.size()-1;j++){
		delayvalue+=delay(temp[j],temp[j+1],fmax,fcurrent);
		
	}
	if(delayvaluebest>delayvalue){
		delayvaluebest=delayvalue;
	}
	delayvalue=0;
}
//cout<<maindelaybest<<" "<<delayvaluebest<<endl;

if(delayvaluebest<maindelaybest){//there is a new best delay
	maindelaybest=delayvaluebest;
}
	//cout<<"maindelaybest:"<<maindelaybest<<"delayvaluebest: "<<delayvaluebest<<endl;
	//cout<<y<<endl;
if(maindelaybest==delayvaluebest){//if not changing until 1000 iterations, stop iterating
	y++;
	//cout<<y<<endl;
}
if(maindelaybest!=delayvaluebest){//reset counter if there is a new best
	y=0;
	//cout<<"--------------------------"<<endl;
}
delayvaluebest=9999999;

//cout<<y<<endl;
//CROSSOVER - problem: 1. redundancy of chromosomes -------------------------------------

for(int i=0;i<chromosomes.size();i+=2){
//cout<<"\n"<<i<<endl;
//cout<<"check";
	bool flag = 0;
	deque<int> ctemp1 = chromosomes[i];
	deque<int> ctemp2 = chromosomes[i+1];
	flag=0;
	if(ctemp1.size()<ctemp2.size()){//first chromosome is the shortest
		int length = ctemp1.size();
		int crossover = rand()%length+1;
		if(length==crossover) crossover--;
		deque<int> kids1;
		deque<int> kids2;
		
		kids1.push_back(ctemp1[0]);
		//
		for(int j=1;j<crossover;j++){//get the genes of first chromosome for kids1
			for(int k=0;k<kids1.size();k++){//check for redundancy of genes/chromosomes
				if((ctemp1[j]==kids1[k])&&(delay(kids1[k],ctemp1[j],fmax, fcurrent)<0)){//if there is redundancy or blocked path, raise flag to ignore the chromosome
					flag =1;
					break;
				}
			}
			if(flag==1) break;
			kids1.push_back(ctemp1[j]);//add constraints
		}
		for(int j=crossover;j<ctemp2.size();j++){//concatenate genes of second chromosome
			for(int k=0;k<kids1.size();k++){//check for redundancy of genes/chromosomes
				if((ctemp2[j]==kids1[k])&&(delay(kids1[k],ctemp2[j],fmax, fcurrent)<0)){//if there is redundancy or blocked path, raise flag to ignore the chromosome
					flag =1;
					break;
				}
			}
			if(flag==1) break;
			kids1.push_back(ctemp2[j]);
		}
		//
		kids2.push_back(ctemp2[0]);
		for(int j=1;j<crossover;j++){//get the genes of first chromosome for kids2
			for(int k=0;k<kids2.size();k++){//check for redundancy of genes/chromosomes
					if((ctemp2[j]==kids2[k])&&(delay(kids2[k],ctemp2[j],fmax, fcurrent)<0)){//if there is redundancy or blocked path, raise flag to ignore the chromosome
						flag =1;
						break;
					}
				}
				if(flag==1) break;
				kids2.push_back(ctemp2[j]);
			}
		for(int j=crossover;j<ctemp1.size();j++){
			for(int k=0;k<kids2.size();k++){//check for redundancy of genes/chromosomes
				if((ctemp1[j]==kids2[k])&&(delay(kids2[k],ctemp1[j],fmax, fcurrent)<0)){//if there is redundancy or blocked path, raise flag to ignore the chromosome
					flag =1;
					break;
				}
			}
			if(flag==1) break;
			kids2.push_back(ctemp1[j]);
		}
		for(int j=0;j<kids1.size();j++){
			for(int k=0;k<kids1.size();k++){
				if(j==k){
					continue;
				}
				if(kids1[j]==kids1[k]){
					flag=1;
				}
			}
		}
		for(int j=0;j<kids2.size();j++){
			for(int k=0;k<kids2.size();k++){
				if(j==k){
					continue;
				}
				if(kids2[j]==kids2[k]){
					flag=1;
				}
			}
		}
		if(flag==1){
			i-=2;
			continue;
		}
		kids.push_back(kids1);
		kids.push_back(kids2);
	}
	
	else{//second chromosome is the shortest
		int length = ctemp2.size();
		int crossover = rand()%length+1;
		if(length==crossover) crossover--;
		deque<int> kids1;
		deque<int> kids2;
		
		//
		kids1.push_back(ctemp2[0]);
		for(int j=1;j<crossover;j++){//get the genes of first chromosome for kids1
			for(int k=0;k<kids1.size();k++){//check for redundancy of genes/chromosomes
				if((ctemp2[j]==kids1[k])&&(delay(kids1[k],ctemp2[j],fmax, fcurrent)<0)){//if there is redundancy or blocked path, raise flag to ignore the chromosome
					flag =1;
					break;
				}
			}
			if(flag==1) break;
			kids1.push_back(ctemp2[j]);//add constraints
		}
		for(int j=crossover;j<ctemp1.size();j++){//concatenate genes of second chromosome
			for(int k=0;k<kids1.size();k++){//check for redundancy of genes/chromosomes
				if((ctemp1[j]==kids1[k])&&(delay(kids1[k],ctemp1[j],fmax, fcurrent)<0)){//if there is redundancy or blocked path, raise flag to ignore the chromosome
					flag =1;
					break;
				}
			}
			if(flag==1) break;
			kids1.push_back(ctemp1[j]);
		}
		//
		kids2.push_back(ctemp1[0]);
		for(int j=1;j<crossover;j++){//get the genes of second chromosome for kids2
			for(int k=0;k<kids2.size();k++){//check for redundancy of genes/chromosomes
					if((ctemp1[j]==kids2[k])&&(delay(kids2[k],ctemp1[j],fmax, fcurrent)<0)){//if there is redundancy or blocked path, raise flag to ignore the chromosome
						flag =1;
						break;
					}
				}
				if(flag==1) break;
				kids2.push_back(ctemp1[j]);
			}
		for(int j=crossover;j<ctemp2.size();j++){
			for(int k=0;k<kids2.size();k++){//check for redundancy of genes/chromosomes
				if((ctemp2[j]==kids2[k])&&(delay(kids2[k],ctemp2[j],fmax, fcurrent)<0)){//if there is redundancy or blocked path, raise flag to ignore the chromosome
					flag =1;
					break;
				}
			}
			if(flag==1) break;
			kids2.push_back(ctemp2[j]);
		}
		for(int j=0;j<kids1.size();j++){
			for(int k=0;k<kids1.size();k++){
				if(j==k){
					continue;
				}
				if(kids1[j]==kids1[k]){
					flag=1;
				}
			}
		}
		for(int j=0;j<kids2.size();j++){
			for(int k=0;k<kids2.size();k++){
				if(j==k){
					continue;
				}
				if(kids2[j]==kids2[k]){
					flag=1;
				}
			}
		}
		if(flag==1){
			i-=2;
			continue;
		}
		kids.push_back(kids1);
		kids.push_back(kids2);
	}
	if(flag==1){
			i-=2;
			continue;
		}

}
//cout<<endl;
//cout<<endl;
/* for(int i=0;i<kids.size();i++){
	deque<int> temp = kids[i];
	for(int j=0;j<temp.size();j++){
		cout<<temp[j]<<" ";
	}
	cout<<endl;
}
 */
//MUTATION------------------------------------------------------------------------
/*  cout<<endl;
cout<<"Mutation!"<<endl;
for(int i=0; i<kids.size(); i++) {
	deque<int> temp=kids[i];
	for(int j=0; j<temp.size(); j++) {
		cout<<temp[j]<<" ";
	}
	cout<<endl;
} */

for(int i=0;i<kids.size();i++){
	deque<int> temp = kids[i];
	double chance = ((double)((double)rand()/rand()));
	int chance1 = chance;
	chance -=chance1;
	bool flag=0;//if flag = 1 , there is a problem, mutation should restart
	if (flag==1) {
		flag=0;
		continue;
	}
	if(temp.size()==2){
		continue;
	}
	//cout<<chance<<endl;
	if(chance<.90){//it will mutate
		int replace = rand()%2;
		if(replace==1){//replace
			int place = rand()%(temp.size()-2) +1;
				for(int j=0;j<temp.size();j++){//if the replacing node is not valid, continue without mutation
		/* 			if(flag==1){
						break;
					} */
						int change = rand()%max.size();//value of node to be change
						for(int k=0; k<temp.size(); k++){//check if random number is not repeating
							if(temp[k]==change){//if repeating, break
								flag=1;
								break;
							}
						}
					if((delay(change,temp[place+1],fmax,fcurrent)<=0)&&(delay(temp[place-1],change,fmax,fcurrent)<=0)&&(flag!=0)){//if delay is invalid
						flag=1;
						break;
					}
					else{//delay is valid
						temp[place] = change;
					}
				}
			for(int j=0;j<temp.size();j++){//check redundancy of mutated chromosome
				for(int k=0;k<temp.size();k++){
					if(j==k){
						continue;
					}
					if(temp[j]==temp[k]){
						flag=1;
					}
				}
			}
			if(flag==1){
				i--;
				continue;
			}
			//cout<<"hue";
			kids[i]=temp;
		}
		
		if(replace==0){//remove
		//cout<<"\n Remove"<<endl;
		
		deque<int> temp1;
		bool flag=0;
		//cout<<delay(temp[place+1],temp[place-1],fmax,fcurrent)<<endl;
		//cout<<temp[place-1]<<" "<<temp[place+1]<<endl;
		//cin>>place;
		for(int l=0;l<temp.size();l++){
			int place = rand()%(temp.size()-2) +1;
			if(flag==1) break;
			if((delay(temp[place+1],temp[place-1],fmax,fcurrent)>0)){//if remove is valid
				flag=1;
				for(int j=0; j<temp.size(); j++) {//copy to new temp for transfer to kids
					if(place==j){
						continue;
					}
						//copy
						temp1.push_back(temp[j]);	
				}
				
			for(int j=0;j<temp1.size();j++){//check redundancy of mutated chromosome
				for(int k=0;k<temp1.size();k++){
					if(j==k){
						continue;
					}
					if(temp1[j]==temp1[k]){
						flag=1;
					}
				}
			}
			if(flag==1){
				i--;
				continue;
			}
			//cout<<"hue";
			kids[i] = temp1;//add to kids
			}
		}
/* 		else{
			i--;
		} */
	
	}
	}
	if(flag==1){
		i--;
		flag=0;
		continue;
	}
	/* else{//if mutation chance is <.95
		continue;
	} */
	

}

//SELECTION----------------------------------------------------------------------

for(int i=0; i<chromosomes.size(); i++) {
	deque<int> temp = chromosomes[i];
	deque<int> temp1 = kids[i];
	double totaldelaykids=0;
	double totaldelaychromosomes=0;
	for(int j=0; j<temp.size()-1; j++) {
		totaldelaychromosomes+=delay(temp[j],temp[j+1],fmax,fcurrent);
	}
	for(int j=0; j<temp1.size()-1; j++) {
		totaldelaykids+=delay(temp1[j],temp1[j+1],fmax,fcurrent);
	}
	//cout<<totaldelaychromosomes<<" "<<totaldelaykids;
	//cout<<"Delay chromosomes: "<<totaldelaychromosomes<<" Delay Kids: "<<totaldelaykids<<endl;
	if(totaldelaykids<=totaldelaychromosomes) {
		chromosomes[i]=kids[i];
		//cout<<"check";
	}
}

kids.clear();
}
cout<<"Number of Iterations: "<<itr<<endl;


/* cout<<endl;
cout<<"New mutation"<<endl;
for(int i=0; i<kids.size(); i++) {
	deque<int> temp=kids[i];
	for(int j=0; j<temp.size(); j++) {
		cout<<temp[j]<<" ";
	}
	cout<<endl;
} */ 

/* for(int i=0;i<kids.size();i++){
	deque <int> temp = kids[i];
	for(int j=0;j<temp.size();j++){
		cout<<temp[j];
	}
	cout<<endl;
}
 */
/* for(int i=0;i<40;i++){
	deque < int > temp = chromosomes[i];
	for(int j=0;j<temp.size();j++){
		cout<<temp[j];
	}
cout<<endl;
} */
double delayvaluebest=9999;
cout<<"Final set of chromosomes:"<<endl;
for(int i=0;i<40;i++){
	deque < int > temp = chromosomes[i];
	double delayvalue=0;
	for(int j=0;j<temp.size()-1;j++){
		delayvalue+=delay(temp[j],temp[j+1],fmax,fcurrent);
		cout<<temp[j]<<" ";
		
	}
	if(delayvaluebest>delayvalue)
		delayvaluebest=delayvalue;
	cout<<temp[temp.size()-1];
	cout<<" --- "<<delayvalue;
	delayvalue=0;
cout<<endl;
}
cout<<"Best Delay Constant: "<<delayvaluebest;
return 0;

}











