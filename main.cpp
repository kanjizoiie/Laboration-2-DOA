#include <cstdlib>
#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <list>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <math.h>

template <typename T>
void sel_sort(T begin, T end) {
    for (T current_outer = begin; current_outer != end; current_outer++) {
        T min = current_outer;
        for (T current_inner = current_outer; current_inner != end; current_inner++) {
            if (*current_inner < *min) {
                min = current_inner;
            }
        }
        if (*min != *current_outer)
            std::swap(*min, *current_outer);
    }
}

template <typename T>
void ins_sort(T begin, T end) {
    // Pick the second element as the current element.
    T current = begin + 1;
    while(current != end) {
        T saved = current;
        // If the saved position is lower than the previous position swap them
        while (saved > begin && saved[-1] > saved[0]) {
            std::swap(saved[-1], saved[0]);
            // Move downwards from the saved value
            saved--;
        }
        current++;
    }
}

template<typename T>
T partition(T begin, T end) {
    T pivot = begin;
    T i = begin - 1;
    T j = end + 1;
    while (true) {
        do {
            i++;
        } while (*i < *pivot);
        do {
            j--;
        } while (*j > *pivot);
        
        if (i >= j) {
            return j;
        }
        std::swap(*j, *i);
    }
}

template<typename T>
T partition_mo3(T begin, T end) {
    int mid_pos = (end - begin) / 2;
    T middle = (begin + mid_pos);
    if (*begin > *middle)
        std::swap(*middle, *begin);
    if (*begin > *(end - 1))
        std::swap(*begin, *(end - 1));
    if (*middle > *(end - 1))
        std::swap(*(end - 1), *middle);
    std::swap(*begin, *middle);
    T pivot = begin;
    T i = begin - 1;
    T j = end + 1;
    
    while (true) {
        do {
            i++;
        } while (*i < *pivot);
        
        do {
            j--;
        } while (*j > *pivot);
        
        if (i >= j) {
            return j;
        }
        
        std::swap(*j, *i);
    }
}

template <typename T>
void quick_sort(T begin, T end) {
    if (begin < end) {
        T p = partition(begin, end);
        quick_sort(begin, p);
        quick_sort(p + 1, end);
    }  
}

template <typename T>
void quick_sort_mo3(T begin, T end) {
    if (begin < end) {
        T p = partition_mo3(begin, end);
        quick_sort_mo3(begin, p);
        quick_sort_mo3(p + 1, end);
    }  
}


struct Timing {
    Timing(int amount, int its, double time, double dev): amount(amount), its(its), time(time), dev(dev) {}
    int amount;
    int its;
    double time;
    double dev;
};

#define amnt_list 5

std::vector<std::vector<Timing>> run_sorts(std::vector<std::vector<int>> lists) {
    std::vector<std::vector<Timing>> result;
    std::vector<Timing> temp;
    for (auto list: lists) {
        std::vector<double> times;
        std::vector<std::vector<int>> current_list(amnt_list, list);
        for (std::vector<int>& inner: current_list) {
            auto start = std::chrono::steady_clock::now();
            quick_sort(inner.begin(), inner.end());
            auto end = std::chrono::steady_clock::now();
            std::chrono::duration<double, std::nano> dur = end - start;
            times.push_back(dur.count());
            std::cout << "Done with " << current_list[0].size() << " elements" << std::endl;
            std::cout << dur.count() << std::endl;
        }
        //Calculate standard deviation
        double mean = std::accumulate(times.begin(), times.end(), 0.0) / (double)times.size();
        double innersum = std::accumulate( times.begin(), times.end(), 0.0, [mean](double a, double b){ return a + ((b - mean) * (b - mean)); });
        double dev = sqrtf((innersum) / ((double)times.size() - 1));
        double total_time = std::accumulate( times.begin(), times.end(), 0.00);
        temp.push_back(Timing(current_list[0].size(), current_list.size(), mean, dev));
    }
    result.push_back(temp);
    temp.clear();

    for (auto list: lists) {
        std::vector<double> times;
        //Create a number of copies of the current list
        std::vector<std::vector<int>> current_list(amnt_list, list);
        //Iterate over the copies of the lists and run a sort that is timed, and the duration calculated and saved
        for (std::vector<int>& inner: current_list) {
            auto start = std::chrono::steady_clock::now();
            quick_sort_mo3(inner.begin(), inner.end());
            auto end = std::chrono::steady_clock::now();
            std::chrono::duration<double, std::nano> dur = end - start;
            times.push_back(dur.count());
            std::cout << "Done with " << current_list[0].size() << " elements" << std::endl;
        }
        //Calculate standard deviation
        double mean = std::accumulate( times.begin(), times.end(), 0.0) / (double)times.size();
        double innersum = std::accumulate( times.begin(), times.end(), 0.0, [mean](double a, double b){ return a + ((b - mean) * (b - mean)); });
        double dev = sqrtf((innersum) / ((double)times.size() - 1));
        double total_time = std::accumulate( times.begin(), times.end(), 0.00);
        temp.push_back(Timing(current_list[0].size(), current_list.size(), mean, dev));
    }
    result.push_back(temp);
    temp.clear();
    
    for (auto list: lists) {
        std::vector<double> times;
        //Create a number of copies of the current list
        std::vector<std::vector<int>> current_list(amnt_list, list);
        //Iterate over the copies of the lists and run a sort that is timed, and the duration calculated and saved
        for (std::vector<int>& inner: current_list) {
            auto start = std::chrono::steady_clock::now();
            sel_sort(inner.begin(), inner.end());
            auto end = std::chrono::steady_clock::now();
            std::chrono::duration<double, std::nano> dur = end - start;
            times.push_back(dur.count());
            std::cout << "Done with " << current_list[0].size() << " elements" << std::endl;
        }
        //Calculate standard deviation
        double mean = std::accumulate( times.begin(), times.end(), 0.0) / (double)times.size();
        double innersum = std::accumulate( times.begin(), times.end(), 0.0, [mean](double a, double b){ return a + ((b - mean) * (b - mean)); });
        double dev = sqrtf((innersum) / ((double)times.size() - 1));
        double total_time = std::accumulate( times.begin(), times.end(), 0.00);
        temp.push_back(Timing(current_list[0].size(), current_list.size(), mean, dev));
    }
    result.push_back(temp);
    temp.clear();
    
    for (auto list: lists) {
        std::vector<double> times;
        //Create a number of copies of the current list
        std::vector<std::vector<int>> current_list(amnt_list, list);
        //Iterate over the copies of the lists and run a sort that is timed, and the duration calculated and saved
        for (std::vector<int>& inner: current_list) {
            auto start = std::chrono::steady_clock::now();
            ins_sort(inner.begin(), inner.end());
            auto end = std::chrono::steady_clock::now();
            std::chrono::duration<double, std::nano> dur = end - start;
            times.push_back(dur.count());
            std::cout << "Done with " << current_list[0].size() << " elements" << std::endl;
        }
        //Calculate standard deviation
        double mean = std::accumulate( times.begin(), times.end(), 0.0) / (double)times.size();
        double innersum = std::accumulate( times.begin(), times.end(), 0.0, [mean](double a, double b){ return a + ((b - mean) * (b - mean)); });
        double dev = sqrtf((innersum) / ((double)times.size() - 1));
        double total_time = std::accumulate( times.begin(), times.end(), 0.00);
        temp.push_back(Timing(current_list[0].size(), current_list.size(), mean, dev));
    }
    result.push_back(temp);
    temp.clear();
    
    for (auto list: lists) {
        std::vector<double> times;
        //Create a number of copies of the current list
        std::vector<std::vector<int>> current_list(amnt_list, list);
        //Iterate over the copies of the lists and run a sort that is timed, and the duration calculated and saved
        for (std::vector<int>& inner: current_list) {
            auto start = std::chrono::steady_clock::now();
            std::sort(inner.begin(), inner.end());
            auto end = std::chrono::steady_clock::now();
            std::chrono::duration<double, std::nano> dur = end - start;
            times.push_back(dur.count());
            std::cout << "Done with " << current_list[0].size() << " elements" << std::endl;            
        }
        //Calculate standard deviation
        double mean = std::accumulate( times.begin(), times.end(), 0.0) / (double)times.size();
        double innersum = std::accumulate( times.begin(), times.end(), 0.0, [mean](double a, double b){ return a + ((b - mean) * (b - mean)); });
        double dev = sqrtf((innersum) / ((double)times.size() - 1));
        double total_time = std::accumulate( times.begin(), times.end(), 0.00);
        temp.push_back(Timing(current_list[0].size(), current_list.size(), mean, dev));
    }
    result.push_back(temp);
    temp.clear();
    
    return result;
}

/*
 * 
 */
int main(int argc, char** argv) {
    //Create the vectors that will be sorted.
    static std::vector<std::vector<int>> lists;
    lists.push_back(std::vector<int>(20000));
    lists.push_back(std::vector<int>(40000));
    lists.push_back(std::vector<int>(60000));
    lists.push_back(std::vector<int>(80000));
    lists.push_back(std::vector<int>(100000));
    lists.push_back(std::vector<int>(200000));
    lists.push_back(std::vector<int>(300000));
    lists.push_back(std::vector<int>(400000));
    lists.push_back(std::vector<int>(500000));
    lists.push_back(std::vector<int>(600000));

    
    
    //Generate random values
    for (auto list: lists) {
        std::generate(list.begin(), list.end(), std::rand);
    }
    
    //Run the sorts on the random values
    std::cout << "RANDOM" << std::endl;
    auto rand = run_sorts(lists);
      
    //Sort the lists
    for (auto& list: lists) {
        std::generate(list.begin(), list.end(), std::rand);
        std::sort(list.begin(), list.end());
        for (auto elem: list) {
            std::cout << elem << std::endl;
        }
    }
    
    //Run the sorts on the sorted lists
    std::cout << "SORTED" << std::endl;
    auto s = run_sorts(lists);
    
    //Reverse the list
    for (auto& list: lists) {
        std::generate(list.begin(), list.end(), std::rand);
        std::sort(list.begin(), list.end());
        std::reverse(list.begin(), list.end());
    }
    
    //Run the sorts on the reversed lists
    std::cout << "REVERSE" << std::endl;
    auto r = run_sorts(lists);
    
    //Create flat lists
    for (auto& list: lists) {
        list = std::vector<int>(list.size(), 5);   
    }
    
    //Run the sorts on the flat lists
    std::cout << "FLAT" << std::endl;
    auto f = run_sorts(lists);
    
    

    std::ofstream random("random.csv");
    std::ofstream sorted("sorted.csv");
    std::ofstream reversed("reversed.csv");
    std::ofstream flat("flat.csv");
    
    random << "length,time,dev,amount" << std::endl;
    for (auto outer: rand) {
        for (auto inner: outer) {
            random  << inner.amount << "," << inner.time << "," << inner.dev << "," << inner.its << std::endl;
        }
        random << std::endl;
    }
    sorted << "length,time,dev,amount" << std::endl;
    for (auto outer: s) {
        for (auto inner: outer) {
            sorted << inner.amount << "," << inner.time << "," << inner.dev << "," <<inner.its << std::endl;
        }
        sorted << std::endl;
    }
    reversed << "length,time,dev,amount" << std::endl;
    for (auto outer: r) {
        for (auto inner: outer) {
            reversed << inner.amount << "," << inner.time << "," << inner.dev << "," << inner.its << std::endl;
        }
        reversed << std::endl;
    }
    flat << "length,time,dev,amount" << std::endl;
    for (auto outer: f) {
        for (auto inner: outer) {
            flat  << inner.amount << "," << inner.time << "," << inner.dev << "," << inner.its << std::endl;
        }
        flat << std::endl;
    }
    
    return 0;
}

