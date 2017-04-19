#include <sdsl/suffix_arrays.hpp>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
using namespace std;

string revComp(string a){

	string b="";
	transform(a.begin(), a.end(),a.begin(), ::toupper);
	reverse(a.begin(),a.end());
	for(int i=0;i < a.size(); i++){
		if(a[i]=='A'){
			b+='T';
		} else if(a[i]=='C'){
			b+='G';
		} else if(a[i]=='G'){
			b+='C';
		} else if(a[i]=='T'){
			b+='A';
		} else {
			b+=a[i];
		}
	}
	return b;

}

int main(int argc, char*argv[]) { 
	/*
		run - ./program reference.fa reads.fa Logfile.txt


	*/
    // sdsl::csa_wt<> fm_index;
    // sdsl::construct_im(fm_index, "mississippi!", 1);
    // cout << "'si' occurs " << sdsl::count(fm_index,"si") << " times.\n";
    // sdsl::store_to_file(fm_index,"fm_index-file.sdsl");
    // ofstream out("fm_index-file.sdsl.html");
    // sdsl::write_structure<sdsl::HTML_FORMAT>(fm_index,out);
    if(argc<4){
        cout << "Not enough arguments, please pass the reference genome, output file and the input file(s)\n";
    }
    vector<string> inputFileNames;
    for(int i=3; i < argc; i++){
        inputFileNames.push_back(argv[i]);
    }
    string alignmentResults=argv[2];
	string delim = "\t";
    string reference = "",a;
    ifstream in;
    in.open(argv[1]);

    // in.open("try.fa");
    while(getline(in,a)){
    	// cout << a << endl;
    	if(a[0]=='>'){
    		reference+='N';
    	}else{
    		a.erase(std::remove(a.begin(), a.end(), '\n'), a.end());
    		transform(a.begin(), a.end(),a.begin(), ::toupper);
    		reference+=a;
    	}
    }
    in.close();
    cout << reference.size() << endl;
    // cout << reference << endl;


    sdsl::csa_wt<> fm_index;
    sdsl::construct_im(fm_index, reference, 1);
    // cout << "'si' occurs " << sdsl::count(fm_index,"si") << " times.\n";
    cout << "Suffix array created\n";

    vector<string> contigs;
    vector<string> segmentId;
    string tempSegmentId="";
    int counterOfContigs = 0;
    string tempContig="";

    // Suffix tree has been created and now we are looking through the input files.


    ofstream alignRes(alignmentResults);

    string tempres = "\%exact\n";

    alignRes << tempres;

    for(int ijk=0; ijk<inputFileNames.size();ijk++){

        cout << "reading - " << inputFileNames[ijk] << endl;
        in.open(inputFileNames[ijk]);
        while(getline(in,a)){
        	// cout << a << endl;
        	if(a[0]=='>'){
        		if(counterOfContigs!=0){
        			contigs.push_back(tempContig);
        		}
                int tempI=1;
                tempSegmentId="";
                while(a[tempI]!=' '){
                    tempSegmentId+=a[tempI];
                    tempI++;
                }
                segmentId.push_back(tempSegmentId);
        		tempContig="";
        		counterOfContigs++;
        	}else{
        		a.erase(std::remove(a.begin(), a.end(), '\n'), a.end());
        		transform(a.begin(), a.end(),a.begin(), ::toupper);
        		tempContig+=a;
        	}
        }
        
        contigs.push_back(tempContig);
        cout << "Read the reads file\n";
        cout << "Total contigs - " << contigs.size() << endl;

        string logFileName = inputFileNames[ijk] + ".exact";
        ofstream out(logFileName);

        int commonContigsCounter=0;

        for(int i=0; i<contigs.size();i++){
        	if(sdsl::count(fm_index,contigs[i]) >= 1 || sdsl::count(fm_index,revComp(contigs[i])) >= 1){
        		commonContigsCounter++;
        		// out << contigs[i] << delim << "1" << endl;
                out << segmentId[i] << delim << "1" << endl;
        	}else{
        		// cout << contigs[i] << endl;
        		// out << contigs[i] << delim << "0" << endl;
                out << segmentId[i] << delim << "0" << endl;
        	}
        }

        float percentage = commonContigsCounter * 100.0 / contigs.size();

        cout << inputFileNames[ijk] << " percentage aligned - " << percentage << endl;

        alignRes << percentage << "\n";
        out.close();
        in.close();

        // cout << revComp("acgtAGAGAGAG") << endl;
    }
    alignRes.close();
}