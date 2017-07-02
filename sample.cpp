#include <iostream>
#include <deque>
#include <cstdlib>
#include <ctime>
#include <cmath>
using namespace std;

deque<double> hue(){
	deque<double> temp;
	for(int i=0;i<100;i++){
		temp.push_back(i);
	}
	return temp;
}

int main(){
	srand(time(NULL));
	for(int i=0;i<10;i++){
		double r = static_cast<double> (rand() / static_cast<double> (RAND_MAX/20));
		cout<<r<<endl;
	}
return 0;
}