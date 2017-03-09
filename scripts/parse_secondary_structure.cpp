#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>


using namespace std;



void findPairs(string structure, vector<int> & pairs)
{
	int i=0, j=0;
 	int len= structure.length();
 	bool found = false;
 	int counter=0;
	while(i < len)
	{
		if( structure[i] =='(' )
		{
			j=i;
			while(!found)
			{				
				if(structure[j] == '(')
					counter++;
				else if(structure[j] == ')')
					counter--;			
				
				if(structure[j] != '.' && counter == 0)
				{
					found = true;
					pairs[i] = j;
					pairs[j] = i;
				}
				j++;
		    }
		    found = false;
		}
        i++;
	}
}

int findNumLoops(string alph, int k, int m) //find the number of hairpin loops or stem branches between k and m
{
	int numLoops=0;
	
	for(int i = k+1;  i < m ; i++)
	{
		if (alph[i] == 'H')
		{
			numLoops++;
			while(alph[i] == 'H')
				i++;		
		}
	}
	return numLoops;
}

/*
takes as input a structure in dot-parantheses format and outputs a string of letters that correspond to elements 
such as hairpin loop, bulge loop etc.
*/
string parse(string structure)
{
	int len = structure.length();
	vector<int> pairs(len, -1);  // create a vector for holding pairs
	findPairs(structure, pairs);// find the pairs 
	
	int i, j, k, m;
	string alph = "";
	string annot = "";
	
	//external bases on the 5' end
	for( i = 0; i < len &&  structure[i] == '.'; i++)
	{
		alph += "E"; //external unpaired region
	}
	
	for( j = i; j < len; j++)
	{
		
        if(structure[j] == '.')
        {
            k=j-1;
        	while(structure[k] == '.')  //search for the nearest parantheses on the left side of this dot
				k--;
			
		    m=j+1;
		    while(m<len && structure[m] == '.') //search for the nearest parantheses on the right side of this dot
		        m++;
		    
		    if(m == len)
		    {				
				annot = "E";
		    }
		    else if(structure[k] == '(' && structure[m] == ')')
            {	
				annot = "H"; //hairpin loop
			}
		    else if(structure[k] == ')' && structure[m] == ')') //bulge or interior
		    {
				if(pairs[m] +1 == pairs[k] )  //m and k's paired bases are consecutive so this "." belongs to a bulge 
					annot = "B"; //hairpin loop
				else  
					annot = "T"; // interior loop
            }
            else if(structure[k] == ')' && structure[m] == '(')
            {
				bool foundEnclosingPair = false;
				//check if there's any i-j pair where i<k and j>m
				for(int pos = 0; pos < k ; pos++)
				{
					if(pairs[pos] > m)
						foundEnclosingPair = true;				
				}
				if(foundEnclosingPair)
					annot = "M"; // multiloop
				else
					annot = "E";
            }
            else if(structure[k] == '(' && structure[m] == '(')  //can be multiloop, interior loop or bulge 
            {
            	if(pairs[m] + 1 == pairs[k])   // m and k's paired bases are consecutive so this "." belongs to a bulge 
					annot = "B";
            	else
            	{
            		annot = "N";  //not known for now, we'll analyze these again below	
            		
                } 
				
            }
            
			while (j < m)  //the neighbor dots belong to the same annotation
			{
		     	alph += annot;
                j++;   
			}
			if(m<len)
			{
		     	if(structure[m] == '(')
         	 		alph += "L"; // now j points to a paired base  
         	 	else if(structure[m] == ')')
         	 		alph +="R";
			}
		}
		else
		{
         	if(structure[j] == '(')
         		alph += "L";
         	else if(structure[j] == ')')
     			alph += "R";
		}
	}
	//analyze multiloops
	for( j = 0; j < len; j++) //get rid of "N"s
	{
		if(alph[j] == 'N')  
		{
			k=j-1;
			while(structure[k] == '.')  //search for the nearest parantheses on the left side of this dot
				k--;
			
			m=j+1;
			while(m<len && structure[m] == '.') //search for the nearest parantheses on the right side of this dot
				m++;
			
			if(findNumLoops(alph, k, m) >= 2 || findNumLoops(alph, m, pairs[k]) >= 2)
				alph[j] = 'M';
			else
				alph[j] = 'T';			
		} 	
	}
	
	return alph;
	
}

/*
Alphabet
L : paired, 5' end  (
R : paired, 3' end  )
H : hairpin loop  
T : internal loop
B : bulge loop
M : multiloop
E : external region
*/

int main(int argc, char * argv[])
{
	ifstream input;
	input.open(argv[1]);
	ofstream out;
	out.open(argv[2]);
	string line;
	string structure;
	
	
	while(!input.eof())
	{
		getline(input, line);
		if(line.find(".") != string::npos)
		{
			istringstream ins;
			ins.str(line);
			ins>>structure;
			string parsedStr = parse(structure);
	//		cout<<"structure "<<structure<<"\t"<<parsedStr<<endl;
			out<<parsedStr<<endl;
		}
	}
	 
	return 0;
}

