#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>
using namespace std;
int main(int argc, char* argv[])
{
	if(argc < 1){
		cout << "Provide inputfile as command line arguments\n";
		cout << "Example - ./a.out input1.txt input2.txt commonKmersFrom1in2.txt commonKmersFrom2in1.txt \n";
		exit(0);
	}
	
	string alignmentResults = argv[2];
	vector<string> inputFileNames,common12,common21;
    for(int i=3; i < argc; i++){
        inputFileNames.push_back(argv[i]);
        string temp = argv[i];
        string abc = ".commonKmers12";
        temp+=abc;
        common12.push_back(temp);
        temp = argv[i];
        abc = ".commonKmers21";
        temp+=abc;
        common21.push_back(temp);
    }
	string inputfile1 = argv[1], inputfile2;
	// cout << "here0\n";
	string outputfile1, outputfile2;
	string km;

	char delim[2]="\t";
	long long int value, count1=0, count2=0, commonKmers=0;

	unordered_map<string, int> kmers;
	unordered_map<string, int> kmersFrom1;
	ifstream input;
	ofstream output;
	// cout << "here1\n";

	input.open(inputfile1);
	if (input.is_open())
	{
		while(input >> km >> value){
			kmers[km]=1;
			count1++;
		}
	}
	input.close();


	ofstream alignRes(alignmentResults);
	string tempRes = "recall";
	tempRes+=delim;
	tempRes+= "precision";
	tempRes+= "\n";
	alignRes << tempRes;



	// cout << "here2\n";
	for(int ijk=0; ijk<inputFileNames.size();ijk++){

		inputfile2 = inputFileNames[ijk];
		outputfile1 = common12[ijk];
		outputfile2 = common21[ijk];
		commonKmers=0;
		count2=0;

		input.open(inputfile2);
		output.open(outputfile1);
		if (input.is_open())
		{
			while(input >> km >> value){
				kmersFrom1[km]=1;
				string c = km + delim;
				if(kmers[km]==1){
					commonKmers++;
					c+='1';
				}else{
					c+='0';
				}
				c+="\n";
				output << c;
				count2++;
			}
		}
		input.close();
		output.close();
		// cout << "here3\n";

		output.open(outputfile2);
		for(auto it=kmers.begin(); it!=kmers.end(); it++){
			string c=it->first + delim;
			if(kmersFrom1[it->first]==1){
				c+='1';
			}else{
				c+='0';
			}
			c+="\n";
			output << c;;
		}
		output.close();
		// cout << "here4\n";

		// cout << count1 << " " << count2 << " " << commonKmers << endl;
		cout << inputfile1 << " total kmers - " << count1 << endl;
		cout << inputfile2 << " total kmers - " << count2 << endl;
		cout <<"common kmers - " << commonKmers << endl;

		double percent = commonKmers*1.0 / count1;
		percent*=100.0;

		cout <<  "common kmers from the total kmers in " <<  inputfile1 << " - "<< percent << " % \n";

		alignRes << percent << "%" << delim;

		percent = commonKmers*1.0 / count2;
		percent*=100.0;

		cout <<  "common kmers from the total kmers in " <<  inputfile2 << " - "<< percent << " % \n";

		alignRes << percent << "%" << endl;

	}

	alignRes.close();





}