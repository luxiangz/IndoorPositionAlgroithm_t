#include <file_operation.h>
#include <string>
#include <iostream>


//************************************
// created:   2017/09/13
// FullName:  CountDataRow
// Access:    public 
// Returns:   int
// Parameter: char * filename

// purpose:	  get the data file ROW
//************************************
int CountDataRow(char *filename)
{
	std::string tmp;
	std::ifstream Readfile;
	Readfile.open(filename, std::ios::in);
	int row=0;//ÐÐ
	if (Readfile.fail())
	{
		return 0;
	} 
	else
	{
		while (getline(Readfile, tmp, '\n'))
		{
			row++;
			//std::cout << tmp << std::endl;
		}
		Readfile.close();
	}
	
	return row;
}
//************************************
// created:   2017/09/13
// FullName:  ReadDataRow
// Access:    public 
// Returns:   std::string
// Parameter: char * filename
// Parameter: int datarow

// purpose:	  
//************************************
std::string ReadDataRow(char *filename, int datarow)
{
	std::string tempcontentrow;
	int i = 0;
	int maxrow=CountDataRow(filename);
	std::ifstream Readfile;
	Readfile.open(filename, std::ios::in);
	
	//int maxrow = CountDataRow(filename);
	
	if (datarow <= 0)
		return "Error 1:ROW Error!\n";
	else if (datarow > maxrow)
		return "Error 2:Out of range!\n";
	else
	{
		while (getline(Readfile, tempcontentrow) && (i < (datarow - 1)))
		{
			//tempcontentrow += tempcontentrow;
			i++;
		}
		Readfile.close();
		return tempcontentrow;
	}
	
}