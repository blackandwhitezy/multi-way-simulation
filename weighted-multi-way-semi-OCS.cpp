/*
input: 
    m: # ground elements
    n: n-way semi-OCS
    N: # experiment repetitions
other variables:
    weight (only using the first n-1 entries)
    status (intermediate vector recording state of each element in a repetition)
    instance (m*m)
    result (m*m*N) (record for each rep, each each round, the state of each element)
    round (realtime intermediate vector recording the configuration in a round)
    input_i, input_j (record selection results over N reps, input_i[x]>0 if i is unselected, -1 if i is selected)
functions:
    setWeight(): Algorithm 1
    setWeight2(): Algorithm 2
    createInstance(): randomly pick n distinct elements for m rounds
    createInstance2(): randomly pick n distinct elements such that no element occurs more than n-1 times
    sampling(): sample according to w(k), w(-1)=0 correspondings to the weight of a selected element
    calcCorr(): statistical expectation in N repetitions
output: 
    corr.txt:
        w(k), instance
        For each repetition, for each round, the atate of each element, and the choice of sampling
        Expectation of correlation for any two elements i,j in each round over N repetitions 
*/

#include <cstdlib> 
#include <stdio.h>
#include <iostream>
#include <vector>
#include <ctime>
#include <bits/stdc++.h> 
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

int m=50;//#ground elements
int n=10;//n-way semi-OCS
int N=1000;//#experiment reps

void setWeight(vector<double> &w,int m,int n){
    w[0]=1;
    for (int i=1; i<m; i++)
        w[i]=w[i-1]*(w[i-1]/(n-1)+1); 
}

void setWeight2(vector<double> &w,int m,int n){
    double a=(double)1/(double)(n-1)*(double)(1-(double)(1)/double(sqrt(n-1)));
    w[0]=1;
    for (int i=1; i<m; i++)
        w[i]=(double)(1)/(double)(1-a*(double)(i)); 
}

void createInstance(vector<vector<int>> &ins, int m, int n){
    vector<int> c(m);//count #appearence of each elements up to each round
    vector<int> v(n);//used to generate n unique random numbers out of m
    int t;
    for (int i=0; i<m; i++){
        for (int j=0;j<m;j++){
            ins[i][j]=-c[j];
        }
        fill(v.begin(),v.end(),-1);
        for (int j=0; j<n; j++){
            t=rand()%m;
            if (find(v.begin(),v.end(),t)!=v.end()){
                j--;
                continue;
            }
            v[j]=t;
            c[t]++;
            ins[i][t]=c[t];
        }
    }
}

void createInstance2(vector<vector<int>> &ins, int m, int n){
    vector<int> c(m);//count #appearence of each elements up to each round
    vector<int> v(n);//used to generate n unique random numbers out of m
    int t;
    for (int i=0; i<floor((n-1)*m/n)-1; i++){
        for (int j=0;j<m;j++){
            ins[i][j]=-c[j];
        }
        fill(v.begin(),v.end(),-1);
        for (int j=0; j<n; j++){
            t=rand()%m;
            if (c[t]>n-2){
                j--;
                continue;
            }
            if (find(v.begin(),v.end(),t)!=v.end()){
                j--;
                continue;
            }
            v[j]=t;
            c[t]++;
            ins[i][t]=c[t];
        }
    }
}

int sampling(vector<double> w, vector<int> tuple){
    vector<double> sums(m);
    double temp=0;
    double target=0;
    int index;
    for (int i=0;i<m;i++){
        if (tuple[i]>-1){
            temp+=w[tuple[i]];
            sums[i]+=temp;
        }
    }
    if (temp==0){
        return -1;
    }
    while(target==(double)(0)){
        target=(double)(rand())/(double)(RAND_MAX);
        target=target*(double)(temp);
    }
    for (int i=0;i<m;i++){
        if (target<=sums[i]){
            index=i;
            break;
        }
    }
    return index;
}

double calcCorr(vector<int> input_i,vector<int> input_j, int N, int ki, int kj){
    double num_and=0;
    double num_i=0;
    double num_j=0;
    double corr;
    for (int t=0;t<N;t++){
        if (input_i[t]>0){
            num_i+=1;
            if (input_j[t]>0){
                num_j+=1;
                num_and+=1;
            }
        }else if(input_j[t]>0){
            num_j+=1;
        }
    }
    if (num_i!=0&&num_j!=0){
        corr=(double)N*num_and/(double)(num_i*num_j);
        cout<<corr<<endl;
    }else{corr=0;} 
    return corr;
}
 



int main(){ 

    /*cout<<"Enter number of ground elements:"<<endl;
    cin>>m;
    cout<<"Enter number of elements per round:"<<endl;
    cin>>n;
    cout<<"Enter number of experiment repetitions:"<<endl;
    cin>>N;*/

    vector<double> weight(m);
    vector<int> status(m);
    vector<vector<int>> instance(m,vector<int>(m));//m*m 
    vector<vector<vector<int>>> result(m,vector<vector<int>>(m,vector<int>(N)));//m*m*N 
    
    srand(time(0));
    ofstream myfile;
    myfile.open ("corr.txt");

    int i,j,l,t,x,ki,kj;
    setWeight2(weight,m,n);
    myfile<<"***WEIGHT***"<<endl;
    for (i=0;i<m;i++){
        myfile<<setw(5)<<i<<setw(20)<<weight[i]<<endl;
    }
    myfile<<endl;

    createInstance2(instance,m,n);
    myfile<<"***INSTANCE***"<<endl;
    for (i=0;i<m;i++){
        for (j=0;j<m;j++){
            myfile<<setw(3)<<j<<":"<<setw(3)<<instance[i][j]<<" ";
        }
        myfile<<endl;
    }
    myfile<<endl;

    int choice;
    vector<int> round(m);
    myfile<<"***EXPERIMENT***"<<endl;
    for (i=0;i<N;i++){
        fill(status.begin(),status.end(),0);
        //myfile<<"The "<<setw(5)<<i<<"th repetition:"<<endl;
        for (j=0;j<m;j++){
            //myfile<<"The "<<setw(5)<<j<<"th round:"<<endl;
            for (x=0;x<m;x++){
                round[x]=instance[j][x]-1;
                if (status[x]<0&&round[x]>-1){
                    round[x]=-1;
                }
                //myfile<<setw(5)<<x<<setw(5)<<round[x]<<endl;
            }
            choice=sampling(weight,round);
            //myfile<<"choice: "<<setw(5)<<choice<<endl;
            if (choice==-1){
                continue;
            }
            status[choice]=-1;
            for (l=0;l<m;l++){
                if (instance[j][l]>0){
                    if (status[l]>-1){
                        status[l]++;
                        //cout<<status[l]<<endl;
                    }
                    //cout<<result[j][l][i]<<endl;
                }
                result[j][l][i]=status[l];
            }
        }
    }

    vector<vector<vector<vector<double>>>> corr(m,vector<vector<vector<double>>>(m,vector<vector<double>>(m,vector<double>(m))));
    vector<int> input_i(N);
    vector<int> input_j(N);
    myfile<<"***CORRELATION***"<<endl;
    for (l=0;l<m;l++){
        for (i=0;i<m;i++){
            if (instance[l][i]>0){
                for (j=i+1;j<m;j++){
                    ki=instance[l][i];
                    kj=abs(instance[l][j]);
                    for (x=0;x<N;x++){
                        input_i[x]=result[l][i][x];
                        input_j[x]=result[l][j][x];
                    }
                    corr[i][j][ki][kj]=calcCorr(input_i,input_j,N,ki,kj);
                }
            }
        }
    }

    double cor;
    for (i=0;i<m;i++){
        for (j=0;j<m;j++){
            for (l=0;l<m;l++){
                for (t=0;t<m;t++){
                    cor=corr[i][j][l][t];
                    if (cor!=(double)0){
                        myfile<<"i: "<<setw(2)<<i<<setw(4)<<"j: "<<setw(2)<<j<<setw(6)<<"ki: "<<setw(2)<<l<<setw(6)<<"kj: "<<setw(2)<<t<<setw(10)<<"corr: "<<setw(20)<<cor<<endl;
                    }
                }
            }
        }
    }

    myfile<<"***POSITIVE CORRELATION***"<<endl;
    for (i=0;i<m;i++){
        for (j=0;j<m;j++){
            for (l=0;l<m;l++){
                for (t=0;t<m;t++){
                    cor=corr[i][j][l][t];
                    if (cor>(double)1.8){
                        myfile<<"i: "<<setw(2)<<i<<setw(4)<<"j: "<<setw(2)<<j<<setw(6)<<"ki: "<<setw(2)<<l<<setw(6)<<"kj: "<<setw(2)<<t<<setw(10)<<"corr: "<<setw(20)<<cor<<endl;
                    }
                }
            }
        }
    }

    myfile.close();
} 