#ifndef _MYFOOD_COMMON_H
#define _MYFOOD_COMMON_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <dirent.h>
#include <vector>
#include <map>
#include <sstream>
#include <iostream>
#include <fstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <malloc.h>

template<typename T>
void output_list(const T &list)
{
    typename T::const_iterator it;
    for(it=list.begin(); it!=list.end(); it++)
    {
        std::cout << *it << ", ";
    }
    std::cout << std::endl << "Size: " << list.size() << std::endl;
}

template<typename T>
void list_to_string(const T&list, const std::string &delimiter, std::string &str)
{
    typename T::const_iterator it;
    std::stringstream ss;
    for(it=list.begin(); it!=list.end(); it++)
    {
        if(it==list.begin())
        {
            ss << *it;
        }
        else
        {
            ss << delimiter << *it;
        }
    }
    str=ss.str();
}

void list_dir(std::string dir_path, std::vector<std::pair<std::string, std::string> > &file_list, std::string sufix_match="")
{
    // return listed file name and full path
    file_list.clear();
    struct dirent* ent = NULL;
    DIR *pDir;
    pDir = opendir((dir_path).c_str());
    if (pDir == NULL) {
        std::cerr << "Cannot list files in directory " << dir_path << std::endl;
        return;
    }
    while (NULL != (ent = readdir(pDir))) {
//        std::cout << int(ent->d_type) << " " << ent->d_name << std::endl;
        std::string dirName(ent->d_name);
        if (strcmp(ent->d_name, ".") == 0 || strcmp(ent->d_name, "..") == 0 || (sufix_match.size()>0 && dirName.find(sufix_match, dirName.size()-sufix_match.size())==std::string::npos)) {
//            std::cout << dirName << " " << dirName.size() << " " << dirName.find(sufix_match, dirName.size()-sufix_match.size()) << " " << dirName.size()-sufix_match.size() << std::endl;
            continue;
        }
        std::string fullDirPath = dir_path + "/" + dirName;
//        std::cout << "File: " << fullDirPath << " " << dirName.find_last_of(sufix_match) << " " << dirName.size() << std::endl;
        file_list.push_back(make_pair(dirName, fullDirPath));
    }
}

void string_split_to_double_list(std::string str, std::string delimiter, std::vector<double> &result)
{
    result.clear();
    size_t last_pos=0, pos=0;
    while((pos=str.find(delimiter, last_pos))!=std::string::npos)
    {
        if(str.substr(pos, delimiter.size())==delimiter)
        {
            result.push_back(atof(str.substr(last_pos, pos-last_pos).c_str()));
            last_pos=pos+delimiter.size();
        }
        else
        {
            last_pos = pos + 1;
        }
    }
    result.push_back(atof(str.substr(last_pos, str.size()-last_pos).c_str()));
}

void string_split_to_int_list(std::string str, std::string delimiter, std::vector<int> &result)
{
    result.clear();
    size_t last_pos=0, pos=0;
    while((pos=str.find(delimiter, last_pos))!=std::string::npos)
    {
        if(str.substr(pos, delimiter.size())==delimiter)
        {
            result.push_back(atoi(str.substr(last_pos, pos-last_pos).c_str()));
            last_pos=pos+delimiter.size();
        }
        else
        {
            last_pos = pos + 1;
        }
    }
    result.push_back(atoi(str.substr(last_pos, str.size()-last_pos).c_str()));
}

void string_split(std::string str, std::string delimiter, std::vector<std::string> &result)
{
    result.clear();
    size_t last_pos=0, pos=0;
    while((pos=str.find(delimiter, last_pos))!=std::string::npos)
    {
        if(str.substr(pos, delimiter.size())==delimiter)
        {
            result.push_back(str.substr(last_pos, pos-last_pos));
            last_pos=pos+delimiter.size();
        }
        else
        {
            last_pos = pos + 1;
        }
    }
    result.push_back(str.substr(last_pos, str.size()-last_pos));
//    std::cout << "Number of sub-strings:" << result.size() << std::endl;
//    output_list<vector<std::string> >(result);
}

template<typename S, typename T>
std::string map_to_json_str(const std::map<S, T> &key_value_map, const std::string &item="")
{
    std::stringstream ss;
    ss << "[\"" << item << "\", {";
    for(typename std::map<S, T>::const_iterator it=key_value_map.begin(); it!=key_value_map.end(); it++)
    {
        if(it==key_value_map.begin())
        {
            ss << "\"" << it->first << "\": " << it->second;   
        }
        else{
            ss << ", \"" << it->first << "\": " << it->second; 
        }
    }
    ss << "}]";
    return ss.str();
}

void display_mallinfo(void)
{
    struct mallinfo mi;
    mi = mallinfo();
    printf("================================================\n");
    printf("Total non-mmapped bytes (arena): %d\n", mi.arena);
    printf("# of free chunks (ordblks):            %d\n", mi.ordblks);
    printf("# of free fastbin blocks (smblks):     %d\n", mi.smblks);
    printf("# of mapped regions (hblks):           %d\n", mi.hblks);
    printf("Bytes in mapped regions (hblkhd):      %d\n", mi.hblkhd);
    printf("Max. total allocated space (usmblks):  %d\n", mi.usmblks);
    printf("Free bytes held in fastbins (fsmblks): %d\n", mi.fsmblks);
    printf("Total allocated space (uordblks):      %d\n", mi.uordblks);
    printf("Total free space (fordblks):           %d\n", mi.fordblks);
    printf("Topmost releasable block (keepcost):   %d\n", mi.keepcost);
    printf("================================================\n");
}

void process_mem_usage(double &vm_usage, double &resident_set)
{
    vm_usage     = 0.0;
    resident_set = 0.0;

    std::ifstream stat_stream("/proc/self/stat");
    std::string pid, comm, state, ppid, pgrp, session, tty_nr;
    std::string tpgid, flags, minflt, cminflt, majflt, cmajflt;
    std::string utime, stime, cutime, cstime, priority, nice;
    std::string O, itrealvalue, starttime;

    // the two fields we want
    //    //
    unsigned long vsize;
    long rss;

    stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
        >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
        >> utime >> stime >> cutime >> cstime >> priority >> nice
        >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest
    
    stat_stream.close();

    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    vm_usage     = vsize / 1024.0;
    resident_set = rss * page_size_kb;    
}

void show_mem_usage()
{
    double vm, rss;
    process_mem_usage(vm, rss);
    printf("[Memory] VM: %f; RSS: %f\n",  vm, rss);
}

#endif
