#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include "kmer_index.h"
using namespace std;


//--------- base_name -----------//__110830__//
void getBaseName(string &in,string &out,char slash,char dot)
{
	int i,j;
	int len=(int)in.length();
	for(i=len-1;i>=0;i--)
	{
		if(in[i]==slash)break;
	}
	i++;
	for(j=len-1;j>=0;j--)
	{
		if(in[j]==dot)break;
	}
	if(j==-1)j=len;
	out=in.substr(i,j-i);
}
void getRootName(string &in,string &out,char slash)
{
	int i;
	int len=(int)in.length();
	for(i=len-1;i>=0;i--)
	{
		if(in[i]==slash)break;
	}
	if(i<=0)out=".";
	else out=in.substr(0,i);
}


//-------- read ATCG -----------//
bool ReadATCG(const char* name, vector<char>& genomes, string &head)
{
    ifstream in(name);
    if(!in.good()) {
        return false;
    }

    //-> skip first header
    string buf;
    if(!getline(in, buf)){
        return false;
    }
    head=buf.substr(1,buf.length()-1);

    //-> read following lines
    while(in.good()){
        char item;
        in>>item;
        if(in.fail()){
            break;
        }
        genomes.push_back(item);
    }
    in.close();

    return true;
}

//---------- genome to signal (DNA) -------------//
void Genomes2SignalSequence(const vector<char>& genomes, 
	vector<int>& index, vector<double>& signals,
	vector<string>& kmer_rec, int scale, 
	int FIVE_or_SIX, int ZSCO_or_NOT, int WANT_TAIL)
{
	long bound;
	if(FIVE_or_SIX==0) //-> 5mer model
	{
		g::Mer2Signal::Genome2Index_5mer(genomes, index);
		bound = genomes.size()-4;//genomes.size()%5;
		signals.assign(bound*scale,0);
		kmer_rec.resize(bound*scale);
		long cur=0;
		for(long i = 0; i < bound; i++){
			//-> get kmer
			string kmer(genomes.begin()+i,genomes.begin()+i+5);
			//-> get signal
			double sigval;
			if(index[i]<0)sigval=0;
			else{
				sigval = g::Mer2Signal::AvgSignalAt_5mer(index[i]);
				if(ZSCO_or_NOT==1) //-> transfer to Zsco
				{
					sigval = (sigval-90.208351)/12.832660;
				}
				else               //-> use original int value
				{
					sigval = (int)(5.7*sigval+14);
				}
			}
			//-> assign
			for(int c = scale; c--;){
				signals[cur]=sigval;
				kmer_rec[cur]=kmer;
				cur++;
			}
		}
	}
	else               //-> 6mer model
	{
		g::Mer2Signal::Genome2Index_6mer(genomes, index);
		bound = genomes.size()-5;//genomes.size()%5;
		signals.assign(bound*scale,0);
		kmer_rec.resize(bound*scale);
		long cur=0;
		for(long i = 0; i < bound; i++){
			//-> get kmer
			string kmer(genomes.begin()+i,genomes.begin()+i+6);
			//-> get signal
			double sigval;
			if(index[i]<0)sigval=0;
			else{
			 	sigval = g::Mer2Signal::AvgSignalAt_6mer(index[i]);
				if(ZSCO_or_NOT==1) //-> transfer to Zsco
				{
					sigval = (sigval-90.208199)/12.868652;
				}
				else               //-> use original int value
				{
					sigval = (int)(5.7*sigval+14);
				}
			}
			//-> assign
			for(int c = scale; c--;){
				signals[cur]=sigval;
				kmer_rec[cur]=kmer;
				cur++;
			}
		}
	}


	//---- tail k_mer ------//
	if(WANT_TAIL==1)
	{
		for(long i = bound; i < genomes.size(); i++){
			for(int c = scale; c--;){
				if(ZSCO_or_NOT==1) //-> transfer to Zsco
				{
					signals.push_back(0);
				}
				else
				{
					signals.push_back(5.7*90.2083+14);
				}
			}
		}
	}

}


//------------ main ------------//
int main(int argc,char **argv)
{
	//---- officail_poremodel ----//__2019_02_06__//
	{
		if(argc<2)
		{
			fprintf(stderr,"./officail_poremodel <genome_fasta> \n");
			exit(-1);
		}
		//input
		string genome_fasta=argv[1];
		int zsco_or_not=0;
		//load
		string name;
		vector<char> genomes;
		ReadATCG(genome_fasta.c_str(), genomes,name);
		//proc
		vector<int> index;
		vector<double> signals;
		vector<string> kmer_rec;
		Genomes2SignalSequence(genomes,index,signals,kmer_rec,1,1,zsco_or_not,1);
		//output
		string outfile="signal_"+name+".kmer";
		FILE *fp=fopen(outfile.c_str(),"wb");
		if(zsco_or_not==1)
		{
			for(long i=0;i<signals.size();i++)fprintf(fp,"%lf\n",signals[i]);
		}
		else
		{
			vector<int> signals_ (signals.begin(), signals.end());
			for(long i=0;i<signals_.size();i++)fprintf(fp,"%d\n",signals_[i]);
		}
		fclose(fp);
		//exit
		exit(0);
	}
}

