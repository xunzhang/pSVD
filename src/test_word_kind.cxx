/* test_word_kind.cpp */

#include <map>
#include <fstream>
#include <vector>
#include <iostream>

int main()
{
	std::map<int, int> word_kind;
	std::vector<int> ac;
	int count = 0;

	std::fstream file;
	file.open("debug5",std::ios::in);
	int aaa;
	while(!file.eof()) {
		file >> aaa;
		ac.push_back(aaa);
	}

	for(int i = 0; i < ac.size(); ++i) {
		if(word_kind.count(ac[i]) == 0) {
			word_kind.insert(std::map<int, int>::value_type(ac[i],1));
			count++;
		}
	}       
	std::cout << "word kind is" << word_kind.size() <<std::endl;
	std::cout << "word kind is" << count <<std::endl;
	file.close();
}

