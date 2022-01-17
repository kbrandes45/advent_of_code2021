
#include <vector>
#include <fstream>
#include "stdio.h"
#include <iostream>
#include <string>
#include <sstream>
#include <set>
#include <numeric>
#include <functional>
#include <bits/stdc++.h>
#include <queue>

void day1() {
    std::ifstream input( "./day1data.txt" );
    if(!input.is_open()){
    std::cerr <<  "file cannot be opened";
    }
    if (!input){
    std::cerr << "errors in file";
    }
    std::vector<int> sea_values;
    std::string line;
    int inc_count = 0;
    while (std::getline(input, line)) {
        int val = std::stoi(line);
        if (!sea_values.empty() && val > sea_values.back()) inc_count++;
        sea_values.push_back(val);
    }
    std::cout<<"Count "<<inc_count<<std::endl;
    int count2 = 0;
    int prev_sum = 0;
    for (std::size_t i = 2; i < sea_values.size(); i++) {
        int summ = sea_values.at(i) + sea_values.at(i-1) + sea_values.at(i-2);
        if (summ > prev_sum && i != 2) count2++;
        prev_sum = summ;
    }
    std::cout<<"Count 2: "<<count2<<std::endl;
}

void forward(const int& change, int* val) {
    *val += change;
}
void down(const int& change, int* val) {
    *val += change;
}
void up(const int& change, int* val) {
    *val -= change;
}
void forwardp2(const int& changex, const int& aim, int* depth, int* x) {
    *x += changex;
    *depth += changex * aim;
}
std::vector<std::string> split_string(std::string& input, const char& splitter) {
    std::stringstream test;
    test<< input;
    std::string segment;
    std::vector<std::string> seglist;

    while(std::getline(test, segment, splitter))
    {
    seglist.push_back(segment);
    }
    return seglist;
}
std::vector<int> split_string_to_int(std::string& input, const char& splitter) {
    std::stringstream test;
    test<< input;
    std::string segment;
    std::vector<int> seglist;

    while(std::getline(test, segment, splitter))
    {
        if (!segment.empty()) seglist.push_back(std::stoi(segment));
    }
    return seglist;
}

void day2() {
    std::ifstream input( "./day2data.txt" );
    if(!input.is_open()){
    std::cerr <<  "file cannot be opened";
    }
    if (!input){
    std::cerr << "errors in file";
    }
    std::vector<std::pair<std::string, int>> directions;
    std::string line;
    int x = 0;
    int depth = 0;
    while (std::getline(input, line)) {
        auto splitseg = split_string(line, ' ');
        if (splitseg.size() != 2) std::cout<<"AHH "<<std::endl;
        directions.push_back(std::make_pair(splitseg.at(0), std::stoi(splitseg.at(1))));
        const auto& dir = directions.back();
        if (dir.first == "forward") forward(dir.second, &x);
        else if (dir.first == "up") up(dir.second, &depth);
        else down(dir.second, &depth);
    }
    std::cout<<"x: "<<x<<" y: "<<depth<<" pdt: "<<x*depth<<std::endl;

    x = 0;
    depth = 0;
    int aim = 0;
    for (const auto& dir : directions) {
        if (dir.first == "forward") forwardp2(dir.second, aim, &depth, &x);
        else if (dir.first == "up") up(dir.second, &aim);
        else down(dir.second, &aim);
    }
    std::cout<<"x: "<<x<<" y: "<<depth<<" aim: "<<aim<<" pdt: "<<x*depth<<std::endl;
}

void day3() {
    std::ifstream input( "./day3data.txt" );
    if(!input.is_open()){
    std::cerr <<  "file cannot be opened";
    }
    if (!input){
    std::cerr << "errors in file";
    }
    std::vector<int> count_of_ones_per_col;
    int binsize = 0;
    int count = 0;
    std::string line = "";
    std::vector<std::string> bitvals;
    while (std::getline(input, line)) {
        if (count_of_ones_per_col.size() == 0) {
            for (std::size_t i = 0; i < line.size(); i++) {
                count_of_ones_per_col.push_back(0);
            }
            binsize = line.size();
        }
        for (std::size_t i = 0; i < line.size(); i++) {
            if (line[i] == '1') count_of_ones_per_col[i]++;
        }
        count++;
        bitvals.push_back(line);
    }
    std::string gamma = "";
    std::string eps = "";
    for (const auto& i : count_of_ones_per_col) {
        if (i > count - i) {
            // majority is 1.
            gamma += "1";
            eps += "0";
        } else {
            gamma += "0";
            eps += "1";
        }
    }
    int gamman = std::stoi(gamma, 0, 2);
    int epsn = std::stoi(eps, 0, 2);
    std::cout<<"G: "<<gamman<<" E: "<<epsn<<" pdt: "<<gamman * epsn<<std::endl;

    std::vector<int> valid_og;
    std::vector<int> valid_co2;
    int s = 0;
    char maj = '1';
    int ogc = 0;
    int co2c = 0;
    while (s < binsize) {

        // keep all with 1s in position s for og count. and the others for co2
        if (s == 0 && valid_og.size() == 0) {
            if (count_of_ones_per_col.at(s) >= count - count_of_ones_per_col.at(s)) {
                maj = '1';
            } else {
                maj = '0';
            }
            for (std::size_t j = 0; j < bitvals.size(); j++) {
                if (bitvals.at(j)[s] == maj) {
                    valid_og.push_back(j);
                } else {
                    valid_co2.push_back(j);
                }
            }
        } else {
            // only go through what each vec cares about.
            std::vector<int> remove_ind;
            std::size_t count_ones = 0;
            for (std::size_t k = 0; k < valid_og.size(); k++) {
                if (bitvals.at(valid_og.at(k)).at(s) == '1') count_ones++;
            }
            maj = (count_ones >= valid_og.size() - count_ones) ? '1' : '0';
            for (std::size_t k = 0; k < valid_og.size(); k++) {
                if (bitvals.at(valid_og.at(k)).at(s) != maj) remove_ind.push_back(k);
                if (valid_og.size() == 2) {
                    std::cout<<" s is "<<s<<" and "<<maj<<" for "<<bitvals.at(valid_og.at(0))<<" or "<<bitvals.at(valid_og.at(1))<<std::endl;
                }
            }
            for (int r = remove_ind.size()-1; r > -1; r--) valid_og.erase(valid_og.begin()+ remove_ind.at(r));
            remove_ind.clear();
            if (valid_og.size() == 1 && ogc == 0) {
                // found it!
                ogc = std::stoi(bitvals.at(valid_og.at(0)), 0, 2);
            }


            count_ones = 0;
            for (std::size_t k = 0; k < valid_co2.size(); k++) {
                if (bitvals.at(valid_co2.at(k)).at(s) == '1') count_ones++;
            }
            maj = (count_ones >= valid_co2.size() - count_ones) ? '1' : '0';
            for (std::size_t k = 0; k < valid_co2.size(); k++) {
                if (bitvals.at(valid_co2.at(k)).at(s) == maj) remove_ind.push_back(k);
            }
            for (int r = remove_ind.size()-1; r > -1; r--) valid_co2.erase(valid_co2.begin()+ remove_ind.at(r));
            remove_ind.clear();
            if (valid_co2.size() == 1 && co2c == 0) {
                // found it!
                co2c = std::stoi(bitvals.at(valid_co2.at(0)), 0, 2);
            }
        }
        s++;
    }
    std::cout<<"co2c "<< co2c<<" og "<<ogc<<" pdt "<<ogc * co2c<<std::endl;

}

int compute_points(const std::vector<std::vector<std::pair<bool, int>>>& board, const int& winner_num) {
    std::size_t sq_size = board.size();
    int summ = 0;
    for (std::size_t i = 0; i < sq_size; i++) {
        for (std::size_t j = 0; j < sq_size; j++) {
            if (!board.at(i).at(j).first) summ += board.at(i).at(j).second;
        }
    }
    return summ * winner_num;
}
bool does_board_win(const std::vector<std::vector<std::pair<bool, int>>>& board) {
    std::size_t sq_size = board.size();
    for (std::size_t i = 0; i < sq_size; i++) {
        std::size_t match_row_count = 0;
        std::size_t match_col_count = 0;
        for (std::size_t j = 0; j < sq_size; j++) {
            if (board.at(i).at(j).first) {
                match_row_count++;}
            if (board.at(j).at(i).first) match_col_count++;
        }
        if (match_row_count == sq_size) return true;
        if (match_col_count == sq_size) return true;
    }
    return false;
}
void update_board(std::vector<std::vector<std::pair<bool, int>>>* board, int drawn_val) {
    std::size_t sq_size = board->size();
    for (std::size_t i = 0; i < sq_size; i++) {
        for (std::size_t j = 0; j < sq_size; j++) {
            if (board->at(i).at(j).second == drawn_val) {
                board->at(i).at(j).first = true;
                return;
            }
        }
    }
}

void day4() {
    std::ifstream input( "./day4data.txt" );
    if(!input.is_open()){
    std::cerr <<  "file cannot be opened";
    }
    if (!input){
    std::cerr << "errors in file";
    }
    std::string line = "";
    std::vector<int> numbers_drawn;
    typedef std::vector<std::vector<std::pair<bool, int>>> Board;
    std::vector<Board> boards;
    Board curr_board;
    while (std::getline(input, line)) {
        if (numbers_drawn.size() == 0) {
            numbers_drawn = split_string_to_int(line, ' ');
        } else if (line.empty()) {
            // If there is a blank line, either working boards is done, or continue.
            if (curr_board.size() == 0) {continue;}
            else {
                boards.push_back(curr_board);
                curr_board.clear();
            }
        } else {
            // make board.
            std::vector<std::pair<bool, int>> row;
            std::vector<int> vals = split_string_to_int(line, ' ');
            for (const auto& v : vals) row.push_back(std::make_pair(false, v));
            curr_board.push_back(row);
        }
    }

    std::set<int> seen;
    for (const auto& num : numbers_drawn) {
        for (std::size_t i = 0; i < boards.size(); i++) {
            update_board(&boards.at(i), num);
            if (does_board_win(boards.at(i))) {
                if (seen.size() == boards.size() - 1 && seen.count(i) == 0) {
                    // LAST ONE TO WIN!
                    int pts = compute_points(boards.at(i), num);
                    std::cout<<"last winners points: "<<pts<<std::endl;
                    return;
                } else {
                    if (seen.size() == 0) {
                        int pts = compute_points(boards.at(i), num);
                        std::cout<<"first winners points "<<pts<<std::endl;
                    }
                    seen.insert(i);
                }
            }
        }
    }

}

void day5() {
    std::ifstream input( "./day5data.txt" );
    if(!input.is_open()){
    std::cerr <<  "file cannot be opened";
    }
    if (!input){
    std::cerr << "errors in file";
    }
    std::string line = "";
    std::vector<int> numbers_drawn;

    // parse strings for segments (split ' ' and take 0/2)
    // find min/max numbers for x and y to make board dimensions.
    // make board: vec<vec<int>> counts
    // loop over lines, find horiz/vert lines, populate board along that dim
    typedef std::pair<int,int> Point;
    std::vector<std::pair<Point, Point>> lines;
    int max_x = 0;
    int max_y = 0;
    while (std::getline(input, line)) {
        auto res = split_string(line, ' ');
        auto pt1 = split_string_to_int(res.at(0), ',');
        auto pt2 = split_string_to_int(res.at(2), ',');
        lines.push_back(std::make_pair(std::make_pair(pt1.at(0), pt1.at(1)), std::make_pair(pt2.at(0), pt2.at(1))));
        if (pt1.at(0) > max_x) max_x = pt1.at(0);
        if (pt2.at(0) > max_x) max_x = pt2.at(0);
        if (pt1.at(1) > max_y) max_y = pt1.at(1);
        if (pt2.at(1) > max_y) max_y = pt2.at(1);
    }
    std::vector<std::vector<int>> counts;
    for (int i = 0; i < max_y+1; i++) {
        std::vector<int> row;
        for (int j = 0; j < max_x+1; j++) {
            row.push_back(0);
        }
        counts.push_back(row);
    }
    bool print = true;
    for (const auto& line : lines) {
        if (line.first.first == line.second.first) {
            // horizontal: fixed x value, moving y value
            int fixed_x = line.first.first;
            int starting_y = (line.first.second <= line.second.second) ? line.first.second : line.second.second;
            int ending_y = (line.first.second >= line.second.second) ? line.first.second : line.second.second;
            for (int i = starting_y; i < ending_y+1; i++) {
                counts.at(i).at(fixed_x)++;
            }
        } else if (line.first.second == line.second.second) {
            // vertical: fixed y value, moving x value.
            int fixed_y = line.first.second;
            int starting_x = (line.first.first <= line.second.first) ? line.first.first : line.second.first;
            int ending_x = (line.first.first >= line.second.first) ? line.first.first : line.second.first;
            for (int i = starting_x; i < ending_x+1; i++) {
                counts.at(fixed_y).at(i)++;
            }
        } else {
            // diagonal line.
            // designate first point as line.first
            // determine sign for change in x and change in y.
            // iterate over magnitude of change.
            int delta_x = std::abs(line.first.first - line.second.first);
            int delta_y = std::abs(line.first.second - line.second.second);
            int sign_x = (line.first.first <= line.second.first) ? 1 : -1;
            int sign_y = (line.first.second <= line.second.second) ? 1 : -1;
            int x = 0;
            int y = 0;
            if (print) {
                std::cout<<"PT1 "<<line.first.first<<" "<<line.first.second<<std::endl;
                std::cout<<"PT2 "<<line.second.first<<" "<<line.second.second<<std::endl;
            }
            while (x < delta_x+1 && y < delta_y+1) {
                int ptx = line.first.first + sign_x * x;
                int pty = line.first.second + sign_y * y;
                if ((x < 10 || delta_x - x < 10 )&& print) std::cout<<"Update "<<ptx<<" "<<pty<<std::endl;
                counts.at(pty).at(ptx)++;
                x++;
                y++;
            }
            if (print) print = false;
        }
    }

    int danger_count = 0;
    for (std::size_t i = 0; i < counts.size(); i++){
        for (std::size_t j = 0; j < counts.at(0).size(); j++) {
            if (counts.at(i).at(j) >= 2) danger_count++;
        }
    }
    std::cout<<"Danger count "<<danger_count<<std::endl;
}

void day6() {
    std::ifstream input( "./day6data.txt" );
    if(!input.is_open()){
    std::cerr <<  "file cannot be opened";
    }
    if (!input){
    std::cerr << "errors in file";
    }
    std::string line = "";
    std::vector<int> latternfish;
    std::vector<long> spawn_day = {0, 0, 0, 0, 0, 0, 0};
    std::vector<long> num_to_wait_on = {0, 0, 0, 0, 0, 0, 0};
    while (std::getline(input, line)) {
        latternfish = split_string_to_int(line, ',');
    }

    for (auto& lat : latternfish) {
        spawn_day.at(lat)++;
    }

    // newbies: track them in num_to_wait_on for their spawn day.

    // spawn every 7 days.
    // new fish start to spawn after 9 days.
    // spawn = add a new fish to the population.
    // spawn a new one at 0 and reset to a counter of 6.
    for (int i = 0; i < 256; i++) {
        int spawn_date = i % 7;
        long new_ones = spawn_day.at(spawn_date) - num_to_wait_on.at(spawn_date);
        // std::cout<<"Spawn day "<<spawn_date<<" with "<<new_ones<<std::endl;
        int new_ones_spawn_date = (spawn_date + 9 )% 7;
        // std::cout<<"New ones spawn date "<<new_ones_spawn_date<<std::endl;
        num_to_wait_on.at(new_ones_spawn_date) = new_ones;
        spawn_day.at(new_ones_spawn_date) += new_ones;
        num_to_wait_on.at(spawn_date) = 0;
    }
    long count = 0;
    for (long& num : spawn_day) count+= num;
    std::cout<<"Count "<<std::fixed<<count<<std::endl;

}
void day7() {
    std::ifstream input( "./day7data.txt" );
    if(!input.is_open()){
    std::cerr <<  "file cannot be opened";
    }
    if (!input){
    std::cerr << "errors in file";
    }
    std::string line = "";
    std::vector<int> horiz;
    while (std::getline(input, line)) {
        horiz = split_string_to_int(line, ',');
    }
    int summ = 0;
    for (const int& h : horiz) summ+= h;
    double avg = std::accumulate(horiz.begin(), horiz.end(), 0.0) /horiz.size();
    int avgi = summ / horiz.size();
    std::cout<<"Here: "<<avg<<" "<<avgi<<std::endl;
    int used_avg = avgi+1;
    if (avg - avgi >= 0.5) used_avg = avgi;
    std::sort(horiz.begin(), horiz.end());
    int half = horiz.size()/2.0;
    int med = horiz.at(half);
    std::cout<<"Avg is "<<used_avg<<" med "<<med<<std::endl;
    double total_change = 0;
    for (const int& h: horiz) {
        int change = 0;
        if (h > used_avg) {
            change = h - used_avg;
        } else if (h < used_avg) {
            change = used_avg - h;
        }
        //d*(d + 1)/2
        total_change += (change * (change+1) / 2.0);
        // for (int i = 0; i < change; i++) {
        //     total_change+= i+1;
        // }
    }
    std::cout<<"Change "<<std::fixed<<total_change<<std::endl;
}

void day8() {
    std::ifstream input( "./day8data.txt" );
    if(!input.is_open()){
    std::cerr <<  "file cannot be opened";
    }
    if (!input){
    std::cerr << "errors in file";
    }
    std::string line = "";
    typedef std::vector<std::string> StringV;
    std::vector<std::pair<StringV, StringV>> parsed;
    while (std::getline(input, line)) {
        std::vector<std::string> s = split_string(line, '|');
        StringV vals = split_string(s.at(1), ' ');
        StringV keys = split_string(s.at(0), ' ');
        parsed.push_back(std::make_pair(keys, vals));
    }
    int count = 0;
    std::set<int> unique_sizes = {2, 3, 4, 7};
    for (const auto& p : parsed) {
        for (const auto& v: p.second) {
            if (unique_sizes.count(v.size()) == 1) count++;
        }
    }
    std::cout<<"Count 1: "<<count<<std::endl;

    // 6 digits: 0, 6, 9. Only 6 will not have values for c/f
    // 5 digits: 2, 3, 5.
    // Determine value of c segment via finding the six digit second without either values from the 2 digit seg
    // Determine value of f by using the other value in the 2 digit seg
    // Determine value of a by using the value in the 3 digit seg (only one 3 digit num).
    // Determine value of d by using 0 and seeing which digit it is missing of remaining two digits in 4.
    // Determine value of b by using 0 and seeing which digit of remaining two digits in 4.
    // 0 will have acf and 6 digits and not one of the two remaining digits of 4 (b,d).
    // 1 will be 2 digits.
    // 2 will be the remaining  5 digit numbers.
    // 3 will be 5 digits and contain all values in the 3 digit number (7)
    // 4 will be 4 digits.
    // 5 will be 5 digits and not contain c.
    // 6 will be af but no c and 6 digits.
    // 7 will be 3 digits.
    // 8 will be 7 digits
    // 9 will be acf and 6 digits and all of the letters in 4.
    int total_sum = 0;
    for (const auto& p : parsed) {
        std::map<int, std::string> num_to_let;
        std::vector<std::string> sixes;
        std::vector<std::string> fives;
        for (const auto& k: p.first) {
            // populate values for unique numbers.
            if (k.size() == 2) num_to_let[1] = k;
            if (k.size() == 3) num_to_let[7] = k;
            if (k.size() == 4) num_to_let[4] = k;
            if (k.size() == 7) num_to_let[8] = k;
            if (k.size() == 5) fives.push_back(k);
            if (k.size() == 6) sixes.push_back(k);
        }
        char letc = 'a';
        for (const auto& s : sixes) {
            // 0 is missing one of the values in 4
            // 6 is missing one of the values of 1.
            // 9 has all values in 4.
            int count_ones = 0;
            for (auto& let : num_to_let[1]) {
                if (s.find(let) != std::string::npos) count_ones++;
                else letc = let;
            }
            if (count_ones == 1) {
                num_to_let[6] = s;
                continue;
            }
            int count_fours = 0;
            for (auto& let : num_to_let[4]) {
                if (s.find(let) != std::string::npos) count_fours++;
            }
            if (count_fours == 3) num_to_let[0] = s;
            if (count_fours == 4) num_to_let[9] = s;
        }
        for (const auto& f : fives) {
            // 5 is missing c.
            // 3 has all values in number 7
            // 2 is the other one.
            int count_sevens = 0;
            for (auto& let : num_to_let[7]) {
                if (f.find(let) != std::string::npos) count_sevens++;
            }
            if (f.find(letc) == std::string::npos) num_to_let[5] = f;
            else if (count_sevens == 3) num_to_let[3] = f;
            else num_to_let[2] = f;
        }

        std::map<std::string, int> let_to_num;
        for (auto& p : num_to_let) {
            std::string v = p.second;
            std::sort(v.begin(), v.end());
            std::cout<<"Let "<< v<<" Num "<<p.first<<std::endl;
            let_to_num[v] = p.first;
        }
        int final_num = 0;
        for (auto& v : p.second) {
            std::cout<<"Hi: "<<v<<std::endl;
        }
        for (int i = 1; i < 5; i++) {
            std::string val_s = p.second.at(i);
            std::sort(val_s.begin(), val_s.end());
            int val = let_to_num[val_s];
            int mult = 1000 / (std::pow(10, i-1));
            std::cout<<"Val "<<val_s<<" with num "<<val<<" mult "<<mult<<std::endl;
            final_num += (mult* val);
        }
        std::cout<<"Final "<<final_num<<std::endl;
        total_sum+=final_num;
    }
    std::cout<<"Total "<<total_sum<<std::endl;

}

std::vector<std::pair<int, int>> get_neighbors(std::vector<std::vector<int>>& parsed, const std::size_t& r, const std::size_t& c) {
    std::vector<std::pair<int, int>> neighbors;
    if (((int)r - 1) >= 0) {
        neighbors.push_back(std::make_pair((int)(r-1), (int)c));
    }
    if(r + 1 < parsed.at(0).size()) {
        neighbors.push_back(std::make_pair((int)(r+1), (int)c));
    }
    if (((int)c - 1) >= 0) {
        neighbors.push_back(std::make_pair((int)r, (int)(c-1)));
    }
    if (c+1 < parsed.size()) {
        neighbors.push_back(std::make_pair((int)r, (int)(c+1)));
    }
    return neighbors;
}

bool check_neighbors(std::vector<std::vector<int>>& parsed, const std::size_t& r, const std::size_t& c) {
    int val = parsed.at(c).at(r);
    if (((int)r - 1) >= 0) {
        int neighv = parsed.at(c).at(r-1);
        if (neighv <= val) return false;
    }
    if(r + 1 < parsed.at(0).size()) {
        int neighv = parsed.at(c).at(r+1);
        if (neighv <= val) return false;
    }
    if (((int)c - 1) >= 0) {
        if (parsed.at(c-1).at(r) <= val) return false;
    }
    if (c+1 < parsed.size()) {
        if (parsed.at(c+1).at(r) <= val) return false;
    }
    return true;
}

void day9() {
    std::ifstream input( "./day9data.txt" );
    if(!input.is_open()){
    std::cerr <<  "file cannot be opened";
    }
    if (!input){
    std::cerr << "errors in file";
    }
    std::string line = "";
    std::vector<std::vector<int>> parsed;
    while (std::getline(input, line)) {
        std::vector<int> s;
        for (const auto& c : line) {
            s.push_back(std::stoi(std::string(1, c)));
        }
        parsed.push_back(s);
    }
    std::vector<std::pair<int,int>> basins;
    int risk_level = 0;
    for (std::size_t r = 0; r < parsed.at(0).size(); r++) {
        for (std::size_t c=0; c < parsed.size(); c++) {
            if (check_neighbors(parsed, r, c)) {
                // lowest!
                risk_level += (parsed.at(c).at(r) + 1);
                basins.push_back(std::make_pair(r,c));
            }
        }
    }
    std::cout<<"Risk level "<<risk_level<<std::endl;
    std::vector<std::vector<bool>> checked;
    for (std::size_t r = 0; r < parsed.size(); r++) {
        std::vector<bool> s;
        for (std::size_t c=0; c < parsed.at(0).size(); c++) {
            s.push_back(false);
        }
        checked.push_back(s);
    }
    std::cout<<" checked dimensions "<<checked.size()<<" "<<checked.at(0).size()<<std::endl;
    std::cout<<" regular dimensions "<<parsed.size()<<" "<<parsed.at(0).size()<<std::endl;
    std::vector<int> basin_sizes;
    // take a basin, expand neighbors, check if in basin, mark expanded
    for (const auto& basin: basins) {
        // make set of unvisited nodes.
        std::set<int> visited;
        std::queue<std::pair<int, int>> next_points;
        next_points.push(basin);
        int basin_count = 0;
        // dfs out from the node until no more points to go.
        while (!next_points.empty()) {
            auto& pt = next_points.front();
            if (!checked.at(pt.second).at(pt.first)) {
                int val = parsed.at(pt.second).at(pt.first);
                checked.at(pt.second).at(pt.first) = true;
                if (val != 9) {
                    // new point for this basin
                    basin_count += 1;
                    auto neighbors = get_neighbors(parsed, pt.first, pt.second);
                    for (const auto& n : neighbors) {
                        if (!checked.at(n.second).at(n.first)) {
                            next_points.push(n);
                        }
                    }
                }
            }
            // remove the pt
            next_points.pop(); // removes first element because queue is FIFO
        }
        basin_sizes.push_back(basin_count);
    }
    int size1 = 0;
    int size2 = 0;
    int size3 = 0;
    for (int bsize : basin_sizes) {
        if (bsize > size3) {
            size1 = size2;
            size2 = size3;
            size3 = bsize;
        } else if (bsize > size2) {
            size1 = size2;
            size2 = bsize;
        } else if (bsize > size1) size1=bsize;
    }
    int size = size1*size2*size3;
    std::cout<<"Size "<<size<<std::endl;
}

void day10() {
    std::ifstream input( "./day10data.txt" );
    if(!input.is_open()){
    std::cerr <<  "file cannot be opened";
    }
    if (!input){
    std::cerr << "errors in file";
    }
    std::string line = "";
    std::vector<std::string> parsed;
    std::set<char> opens = {'{', '[', '<', '('};
    std::set<char> closed = {'}', ']', '>', ')'};
    std::map<char,char> matches = {{'{','}'}, {'[',']'}, {'<','>'}, {'(',')'}};
    std::map<char, int> points_values = {{'}',1197}, {']',57},{')',3}, {'>',25137}};
    std::map<char, int> noncorrupt_point_values = {{'}',3}, {']',2},{')',1}, {'>',4}};
    int points = 0;
    std::vector<long> noncorrupt_scores;
    while (std::getline(input, line)) {
        bool is_corrupt= false;
        std::vector<char> open_characters;
        for (const auto& c : line) {
            if (opens.count(c) == 1) open_characters.push_back(c);
            else {
                // get the last open one and make sure they match.
                if (matches[open_characters.back()] == c) {
                    open_characters.pop_back();
                }
                else {
                    // illegal
                    points += points_values[c];
                    is_corrupt = true;
                    break;
                }
            }
        }
        if (!is_corrupt && open_characters.size() > 0) {
            // iterate in reverse order through remaining chars
            long non_corrupt_line_score = 0;
            for (std::size_t i = open_characters.size(); i > 0; i--) {
                int index = i - 1;
                non_corrupt_line_score = 5 * non_corrupt_line_score + noncorrupt_point_values[matches[open_characters.at(index)]];
            }
            noncorrupt_scores.push_back(non_corrupt_line_score);
            std::cout<<"line score "<<non_corrupt_line_score<<std::endl;
        }
    }
    std::cout<<"Illegal points: "<<points<<std::endl;
    std::sort(noncorrupt_scores.begin(), noncorrupt_scores.end());
    int index = (noncorrupt_scores.size() - 1) / 2; // works because it should be zero indexed!
    std::cout<<"Non corrupt score "<<noncorrupt_scores.at(index)<<std::endl;
}

std::vector<std::pair<int,int>> flash_and_return_next(std::vector<std::vector<int>>& parsed, int i, int j) {
    std::vector<std::pair<int, int>> neighbors = {{i+1, j}, {i+1, j+1}, {i+1, j-1}, {i-1, j}, {i-1, j-1}, {i-1, j+1}, {i, j-1}, {i, j+1}};
    std::vector<std::pair<int, int>> next_flashes;
    for (auto& n : neighbors) {
        if (n.first < 0 || n.first >= (int)parsed.size()) continue;
        if (n.second < 0 || n.second >= (int)parsed.at(0).size()) continue;
        parsed.at(n.first).at(n.second)++;
        if (parsed.at(n.first).at(n.second) > 9) next_flashes.push_back(n);
    }
    return next_flashes;
}

void day11() {
    std::ifstream input( "./day11data.txt" );
    if(!input.is_open()){
    std::cerr <<  "file cannot be opened";
    }
    if (!input){
    std::cerr << "errors in file";
    }
    std::string line = "";
    std::vector<std::vector<int>> parsed;
    std::vector<std::vector<bool>> flashed_this_round;
    while (std::getline(input, line)) {
        std::vector<int> s;
        std::vector<bool> b;
        for (const auto& c : line) {
            s.push_back(std::stoi(std::string(1, c)));
            b.push_back(false);
        }
        flashed_this_round.push_back(b);
        parsed.push_back(s);
    }
    int total_flashes = 0;
    // for (int i = 0; i < 100; i++) {
    bool found = false;
    int iter = 0;
    int octo_size = flashed_this_round.size() * flashed_this_round.at(0).size();
    while (!found) {
        // increase everything by 1. Identifying any that are greater than 9
        std::queue<std::pair<int, int>> to_flash;
        for (std::size_t i = 0; i < parsed.size(); i++) {
            for (std::size_t j = 0; j < parsed.at(0).size(); j++) {
                parsed.at(i).at(j)++;
                if (parsed.at(i).at(j) > 9) {
                    // flasher.
                    to_flash.push(std::make_pair(i, j));
                }
            }
        }
        // after the initial round, keep flashing until no more left.
        int fc = 0;
        while (!to_flash.empty()) {
            fc++;
            auto& pt = to_flash.front();
            if (flashed_this_round.at(pt.first).at(pt.second)) {
                to_flash.pop();
                continue;
            }
            // update its neighbors and return any one that needs to now flash.
            auto next_flashes = flash_and_return_next(parsed, pt.first, pt.second);
            for (auto& nextp : next_flashes) {
                // check if it already flashed. only add it if its new.
                if (!flashed_this_round.at(nextp.first).at(nextp.second)) {
                    // std::cout<<"next "<<nextp.first<<" "<<nextp.second<<std::endl;
                    to_flash.push(nextp);
                }
            }
            // mark the square as flashed.
            flashed_this_round.at(pt.first).at(pt.second) = true;
            // remove point from list.
            to_flash.pop();
        }
        std::cout<<"init "<<fc<<std::endl;

        // check how many flashed this round and reset the flashed ones to 0.
        int flash_count = 0;
        std::cout<<"ROUND "<<iter<<" ==============="<<std::endl;
        for (std::size_t i = 0; i < parsed.size(); i++) {
            std::string row = "";
            for (std::size_t j = 0; j < parsed.at(0).size(); j++) {
                if (flashed_this_round.at(i).at(j)) {
                    flash_count++;
                    parsed.at(i).at(j) = 0;
                }
                row+=std::to_string(parsed.at(i).at(j));
                // reset things!
                flashed_this_round.at(i).at(j) = false;
            }
            std::cout<<row<<std::endl;
        }
        if (flash_count == octo_size) {
            found = true;
            break;
        }
        total_flashes+= flash_count;
        iter++;
    }
    std::cout<<"Days "<<iter<<std::endl;
    std::cout<<"Total flashes: "<<total_flashes<<std::endl;
}

std::vector<std::string> copy_path(const std::vector<std::string>& path) {
    std::vector<std::string> out_path;
    for (const std::string& node : path) out_path.push_back(node);
    return out_path;
}

bool is_node_in_path(const std::vector<std::string>& path, const std::string& in_node) {
    for (const std::string& node : path) {
        if (node == in_node) return true;
    }
    return false;
}

bool isUpper(const std::string& s) {
    return std::all_of(s.begin(), s.end(), [](unsigned char c){ return std::isupper(c); });
}

bool has_lower_case_cycle_already(const std::vector<std::string>& path) {
    std::map<std::string, int> node_to_count;
    for (const std::string& node : path) {
        if (node_to_count.count(node) > 0) node_to_count.at(node)++;
        else node_to_count[node] = 1;
    }
    for (const auto& nodep : node_to_count) {
        if (!isUpper(nodep.first) && nodep.second > 1) return true;
    }
    return false;
}

void day12() {
    std::ifstream input( "./day12data.txt" );
    if(!input.is_open()){
    std::cerr <<  "file cannot be opened";
    }
    if (!input){
    std::cerr << "errors in file";
    }
    std::string line = "";
    std::map<std::string, std::vector<std::string>> node_to_neighbors;
    while (std::getline(input, line)) {
        std::vector<std::string> path = split_string(line, '-');
        if (node_to_neighbors.count(path.at(0)) > 0) {
            node_to_neighbors.at(path.at(0)).push_back(path.at(1));
        } else {
            node_to_neighbors[path.at(0)] = {path.at(1)};
        }
        // add the reverse direction edge also.
        if (node_to_neighbors.count(path.at(1))> 0) {
            node_to_neighbors.at(path.at(1)).push_back(path.at(0));
        } else {
            node_to_neighbors[path.at(1)] = {path.at(0)};
        }
    }

    // bfs and build out possible paths
    std::queue<std::vector<std::string>> queue;
    std::vector<std::vector<std::string>> paths;
    queue.push({"start"});
    while (queue.size() > 0) {
        // Pop a partial path.
        // take last node of the partial path and look at all of its neighbors.
        // if neighbor is lowercase and is already in the path, then skip that neighbor.
        // if neighbor is lowercase and is not in path, add it to the path.
        // if neighbor is uppercase, add it to the path.
        std::vector<std::string>& partial_path = queue.front();
        std::string& last_node = partial_path.back();
        for (const std::string& neighbor : node_to_neighbors[last_node]) {
            if (neighbor == "end") {
                // found a path that reaches the goal!
                auto new_path = copy_path(partial_path);
                new_path.push_back("end");
                paths.push_back(new_path);
            } else if (isUpper(neighbor)) {
                // we can re-explore it! yay!
                auto new_path = copy_path(partial_path);
                new_path.push_back(neighbor);
                queue.push(new_path);
            } else {
                // it is lower case.
                if (!is_node_in_path(partial_path, neighbor)) {
                    // fully new path, first time exploring this lower case node.
                    auto new_path = copy_path(partial_path);
                    new_path.push_back(neighbor);
                    queue.push(new_path);
                } else {
                    // we have seen this node before! if there are no other single path cycles, then
                    // add in a cycle here.
                    if (!has_lower_case_cycle_already(partial_path) && neighbor != "start" && neighbor != "end") {
                        auto new_path = copy_path(partial_path);
                        new_path.push_back(neighbor);
                        queue.push(new_path);
                    }
                }
            }
        }
        // erase the partial path we just explored.
        queue.pop();
    }

    std::cout<<"Paths: "<<paths.size()<<std::endl;
    // for (const auto& path : paths) {
    //     std::string p = "";
    //     for (const auto& node : path) {
    //         p+= node;
    //         p+= "-> ";
    //     }
    //     std::cout<<"Path: "<<p<<std::endl;
    // }
}

struct Point {
    int x;
    int y;
    int round_it_moved_on = -1;
    Point(int xpos, int ypos) {
        x = xpos;
        y = ypos;
    }
};

bool is_duplicated(const int& newx, const int& newy, const std::vector<Point>& points) {
    for (const Point& p : points) {
        if (p.x == newx && p.y == newy) return true;
    }
    return false;
}

int apply_x_dir(const int& val, Point* point) {
    if (point->x > val) {
        return 2*val - point->x;
    }
    return -1;
}

int apply_y_dir(const int& val, Point* point) {
    if (point->y > val) {
        return 2*val - point->y;
    }
    return -1;
}

void draw(const std::vector<Point>& points) {
    std::cout<<"DRAWING TIME!"<<std::endl;
    // draw the valid points remaining.
    // find min/max to get the range.
    int max_x = 0;
    int max_y = 0;
    std::map<int, std::vector<int>> row_to_cols_drawn;
    for (const auto& pt : points) {
        if (pt.x > max_x) max_x = pt.x;
        if (pt.y > max_y) max_y = pt.y;
        if (row_to_cols_drawn.count(pt.y) > 0) row_to_cols_drawn.at(pt.y).push_back(pt.x);
        else row_to_cols_drawn[pt.y]= {pt.x};
    }

    max_x += 2;
    max_y += 2;

    for (int y = 0; y < max_y; y++) {
        std::string y_row = "";
        std::vector<int> xs;
        if (row_to_cols_drawn.count(y) == 1) {
            xs = row_to_cols_drawn.at(y);
        }
        for (int x = 0; x < max_x; x++) {
            if (std::find(xs.begin(), xs.end(), x) != xs.end()) y_row+="#";
            else y_row+=".";
        }
        std::cout<<y_row<<std::endl;
    }
}

void day13() {
    std::ifstream input( "./day13data.txt" );
    if(!input.is_open()){
    std::cerr <<  "file cannot be opened";
    }
    if (!input){
    std::cerr << "errors in file";
    }
    std::string line = "";

    std::vector<Point> points;
    std::vector<std::pair<char, int>> directions;
    bool is_point_data = true;
    while (std::getline(input, line)) {
        if (line == "") {
            is_point_data = false;
            continue;
        }
        if (is_point_data) {
            std::vector<int> pointdata = split_string_to_int(line, ',');
            points.push_back(Point(pointdata.at(0), pointdata.at(1)));
        } else {
            std::vector<std::string> dirs = split_string(line, '=');
            int val = std::stoi(dirs.at(1));
            if (dirs.at(0).find('x') != std::string::npos) directions.push_back(std::make_pair('x', val));
            else directions.push_back(std::make_pair('y', val));
        }
    }

    for (std::size_t count = 0; count < directions.size(); count++) {
        auto& dir = directions.at(count);
        std::vector<int> points_to_erase;
        for (std::size_t i = 0; i < points.size(); i++) {
            if (dir.first == 'x') {
                int new_x = apply_x_dir(dir.second, &points.at(i));
                if (new_x != -1) {
                    // point has moved.
                    if (is_duplicated(new_x, points.at(i).y, points)) {
                        points_to_erase.push_back(i);
                    }
                    points.at(i).x = new_x;
                }
            } else {
                int new_y = apply_y_dir(dir.second, &points.at(i));
                if (new_y != -1) {
                    if (is_duplicated(points.at(i).x, new_y, points)) {
                        points_to_erase.push_back(i);
                    }
                    points.at(i).y = new_y;
                }
            }

        }
        // remove duplicates.
        for (std::size_t i = points_to_erase.size(); i > 0; i--) {
            points.erase(points.begin() + points_to_erase.at(i - 1));
        }
        std::cout<<"Visible points after "<<count<<" instructions: "<<points.size()<<std::endl;
    }
    draw(points);

}

void day14() {
    std::ifstream input( "./day14data.txt" );
    if(!input.is_open()){
    std::cerr <<  "file cannot be opened";
    }
    if (!input){
    std::cerr << "errors in file";
    }
    std::string line = "";
    std::vector<char> templ;
    std::map<std::string, char> val_to_insert;
    while (std::getline(input, line)) {
        if (line == "") continue;
        std::vector<std::string> results = split_string(line, ' ');
        if (results.size() == 1) {
            // found the template.
            for (std::size_t i = 0; i < results.at(0).size(); i++) {
                templ.push_back(results.at(0).at(i));
            }
        } else {
            std::cout<<"Val to insert: "<<results.at(0)<<std::endl;
            val_to_insert[results.at(0)] = results.at(2).at(0);
        }
    }
    bool part1 = false;
    if (part1) {
        for (int k = 0; k < 10; k++) {
            std::map<int, char> inserts;
            for (std::size_t j = 0; j < templ.size() - 1; j++) {
                std::string one(1, templ.at(j));
                std::string two(1, templ.at(j+1));
                std::string pattern = one+two;
                if (val_to_insert.count(pattern) > 0) {
                    inserts[j] = val_to_insert.at(pattern);
                }
            }
            std::vector<char> updated_templ;
            std::string total_S = "";
            for (std::size_t i = 0; i < templ.size(); i++) {
                updated_templ.push_back(templ.at(i));
                total_S+= templ.at(i);
                if (inserts.count(i) > 0) {
                    updated_templ.push_back(inserts[i]);
                    total_S+= inserts.at(i);
                }
            }
            templ = updated_templ;
            // std::cout<<"Templ at "<<k<<" is "<<total_S<<std::endl;
        }
        std::map<char, long> char_count;
        for (const char& c : templ) {
            if (char_count.count(c) > 0) char_count.at(c)++;
            else char_count[c] = 1;
        }
        long minv = 1000000000000000;
        long maxv = 0;
        for (const auto& p : char_count) {
            if (p.second > maxv) {maxv = p.second; std::cout<<p.first<<" set max at "<<maxv<<std::endl;}
            if (p.second < minv) {minv = p.second; std::cout<<"AND "<<p.first<<" set min at "<<minv<<std::endl;}
        }
        std::cout<<"Diff "<<maxv - minv<<" for "<<maxv<<" "<<minv<<std::endl;
    }

    std::map<std::string, int> upd_partners;
    std::map<std::string, int> adjustments;
    std::string first_string = "";
    std::string last_string = "";
    for (std::size_t i = 0; i < templ.size() - 1; i++) {
        std::string one(1, templ.at(i));
        std::string two(1, templ.at(i+1));
        if (i == 0){
            first_string = one+two;
        } else if (i == templ.size() - 2) {
            last_string = one+two;
        } else {
            upd_partners[one+two] = 1;
            adjustments[one+two] = 0;
        }
    }


    for (int k = 0; k < 2; k++) {
        for (auto& partner_pair : upd_partners) {
            if (partner_pair.second == 0) continue;
            std::cout<<"or deep within the beast? "<<partner_pair.first<<" "<<partner_pair.second<<std::endl;
            std::string val(1, val_to_insert.at(partner_pair.first));
            std::string one(1, partner_pair.first.at(0));
            std::string two(1, partner_pair.first.at(1));
            std::string one_v = one+val;
            std::string two_v = val+two;
            std::cout<<"Creates "<<one_v<<" and "<<two_v<<std::endl;
            // remove it from partners map.
            adjustments[partner_pair.first] -= partner_pair.second;
            std::cout<<"Adjustment for "<<partner_pair.first<<" is "<<adjustments[partner_pair.first]<<std::endl;
            // add in the two new patterns.
            if (adjustments.count(one_v) > 0) adjustments[one_v] += partner_pair.second;
            else adjustments[one_v] = 1;
            std::cout<<"And then adj for "<<one_v<<" is now "<<adjustments[one_v]<<std::endl;
            std::cout<<"prior for "<<two_v<<" is "<<adjustments[two_v]<<std::endl;
            if (adjustments.count(two_v) > 0) adjustments[two_v] += partner_pair.second;
            else adjustments[two_v] = 1;
            std::cout<<"And then adj for "<<two_v<<" is now "<<adjustments[two_v]<<std::endl;
        }
        std::cout<<"this?"<<std::endl;
        // propegate and update first/last strings also.
        std::string val(1, val_to_insert.at(first_string));
        std::string one(1, first_string.at(0));
        std::string two(1, first_string.at(1));
        if (upd_partners.count(first_string) > 0) upd_partners[first_string] -=1;
        else upd_partners[first_string] = 0;
        first_string = one+val;
        std::cout<<" first one "<<val+two<<" with prior adj "<<adjustments[val+two]<<std::endl;
        if (upd_partners.count(val+two) > 0) upd_partners[val+two] += 1;
        else upd_partners[val+two] = 1;
        std::cout<<"here? "<<first_string<<std::endl;
        std::cout<<"last string?"<<last_string<<std::endl;
        std::string val1(1, val_to_insert.at(last_string));
        std::string one1(1, last_string.at(0));
        std::string two1(1, last_string.at(1));
        std::cout<<"and "<<one1+val1<<std::endl;
        if (upd_partners.count(one1+val1) > 0) upd_partners[one1+val1] += 1;
        else upd_partners[one1+val1] = 1;
        if (upd_partners.count(last_string) > 0) upd_partners[last_string] -= 1;
        else upd_partners[last_string] = 0;
        std::cout<<"after subtraction "<<adjustments[last_string]<<std::endl;
        last_string = val1 + two1;
        std::cout<<"UPD last "<<last_string<<std::endl;

        for (auto& adj : adjustments) {
            // apply each adjustment.
            std::cout<<"Adjustment "<<adj.first<<" of "<<adj.second<<std::endl;
            if (upd_partners.count(adj.first) > 0) upd_partners[adj.first] += adj.second;
            else upd_partners[adj.first] = adj.second;
            adjustments[adj.first] = 0;
        }
        std::cout<<"end of fun====================="<<std::endl;
    }

    std::map<char, long> counts;
    for (const auto& partner : upd_partners) {
        std::cout<<"Partner "<<partner.first<<" with count "<<partner.second<<std::endl;
        if (counts.count(partner.first.at(0)) > 0) counts.at(partner.first.at(0)) += partner.second;
        else counts[partner.first.at(0)] = partner.second;
        if (counts.count(partner.first.at(1)) > 0) counts.at(partner.first.at(1)) += partner.second;
        else counts[partner.first.at(1)] = partner.second;
    }

    long minv = 1000000000000000;
    long maxv = 0;
    for (const auto& p : counts) {
        long v = p.second / 2;
        if (p.first == first_string.at(0)) v++;
        // if (p.first == first_string.at(1)) v++;
        // if (p.first == last_string.at(0)) v++;
        if (p.first == last_string.at(1)) v++;
        if (v > maxv) {maxv = v;}// std::cout<<p.first<<" set max at "<<maxv<<std::endl;}
        if (v < minv) {minv = v;}// std::cout<<"AND "<<p.first<<" set min at "<<minv<<std::endl;}
        std::cout<<"Char "<<p.first<<" has "<<v<<std::endl;

    }
    std::cout<<"Diff "<<maxv - minv<<" for "<<maxv<<" "<<minv<<std::endl;
}

double heuristic(std::pair<int, int>& goal, std::pair<int, int> current) {
    // return manhattan distance.
    return std::abs(goal.first - current.first) + std::abs(goal.second - current.second);
}

double cost(std::pair<int, int>& next_node, const std::vector<std::vector<int>>& risk_levels) {
    return risk_levels.at(next_node.first).at(next_node.second);
}

struct RiskPoint {
    int x;
    int y;
    double risk;

    RiskPoint* parent = nullptr;
    double cost_for_path_to_node = std::numeric_limits<double>::max();
    double f_score = std::numeric_limits<double>::max();

    RiskPoint(int in_x, int in_y, double in_risk) {
        x=in_x;
        y=in_y;
        risk = in_risk;
    }

    std::vector<std::pair<int, int>> get_neighbors(const int& max_x, const int& max_y) {
        std::vector<std::pair<int, int>> ret;
        if (x-1 >= 0) ret.push_back(std::make_pair(x-1, y));
        if (x+1 < max_x) ret.push_back(std::make_pair(x+1, y));
        if (y-1 >= 0) ret.push_back(std::make_pair(x, y-1));
        if (y+1 < max_y) ret.push_back(std::make_pair(x,y+1));
        return ret;
    }
};

bool found_point_at(std::vector<std::pair<int, int>>& points, const int& x, const int& y, int* index) {
    for (std::size_t i = 0; i < points.size(); i++){
        if (points.at(i).first == x && points.at(i).second == y) {
            if (index) *index = (int) i;
            return true;
        }
    }
    return false;
}

void astar(std::vector<std::vector<RiskPoint>>& risk_levels) {
    int max_x = risk_levels.size();
    int max_y = risk_levels.at(0).size();
    std::pair<int, int> goal = std::make_pair(max_x-1, max_y-1);

    // track path and its accumulated f-score. Expand the next lowest fsocre each time.
    std::vector<std::pair<int, int>> next_path;
    std::vector<std::pair<int, int>> closed;
    RiskPoint& start = risk_levels.at(0).at(0);
    start.parent = nullptr;
    start.cost_for_path_to_node = 0.0;// start.risk;
    start.f_score = 0.0;
    next_path.push_back(std::make_pair(0,0));
    while (!next_path.empty()) {
        // get lowest f_Score point and remove it.
        int low_index = -1;
        double lowest_f = std::numeric_limits<double>::max();
        for (std::size_t i =0; i < next_path.size(); i++) {
            std::pair<int, int>& pt = next_path.at(i);
            RiskPoint& rp = risk_levels.at(pt.first).at(pt.second);
            if (rp.f_score < lowest_f) {
                lowest_f = rp.f_score;
                low_index = i;
            }
        }
        std::pair<int, int>& nextt = next_path.at(low_index);
        RiskPoint& next_pt = risk_levels.at(nextt.first).at(nextt.second);
        if (nextt.first == goal.first && nextt.second == goal.second) {
            // at the goal!
            std::cout<<"GOAL! "<<next_pt.cost_for_path_to_node<<std::endl;
            return;
        }

        // explore neighbors of the current node.
        auto neighbors = next_pt.get_neighbors(max_x, max_y);
        for (std::pair<int, int>& n_pt : neighbors) {
            RiskPoint& n = risk_levels.at(n_pt.first).at(n_pt.second);
            // get cost from node-> neighbor == risk level at neighbor.
            double node_to_neighbor_cost = n.risk;
            // compute cost of path up to neighbor through this node.
            double cost_to_neighbor_thru_node = next_pt.cost_for_path_to_node + node_to_neighbor_cost;
            // check if this cost is lower than the neighbors exisitng cost.
            if (cost_to_neighbor_thru_node < n.cost_for_path_to_node) {
                // better (or just first) path found!
                n.cost_for_path_to_node = cost_to_neighbor_thru_node;
                n.parent = &next_pt;

                // add it in to the paths to explore if it isn't already there.
                if (!found_point_at(next_path, n.x, n.y, nullptr)) {
                    next_path.push_back(std::make_pair(n.x, n.y));
                }

                // update the fscore for the path.
                n.f_score = cost_to_neighbor_thru_node + heuristic(goal, std::make_pair(n.x, n.y));

                // remove it from closed if it is there.
                int index_to_remove = -1;
                bool found = found_point_at(closed, n.x, n.y, &index_to_remove);
                if (found && index_to_remove != -1) {
                    closed.erase(closed.begin() + index_to_remove);
                }
            }
        }

        // close out this node.
        closed.push_back(nextt);
        next_path.erase(next_path.begin() + low_index);
    }
    std::cout<<"No path?"<<std::endl;
}

void day15() {
    std::ifstream input( "./day15data.txt" );
    if(!input.is_open()){
    std::cerr <<  "file cannot be opened";
    }
    if (!input){
    std::cerr << "errors in file";
    }
    std::string line = "";
    std::vector<std::vector<RiskPoint>> risk_levels;
    while (std::getline(input, line)) {
        // std::vector<int> results = split_string_to_int(line, ',');
        std::vector<RiskPoint> row_levels;
        for (std::size_t i = 0; i < line.size(); i++) {
                std::stringstream test;
            std::string c(1, line[i]);
            row_levels.push_back(RiskPoint((int)risk_levels.size(), i, std::stoi(c)));
        }
        risk_levels.push_back(row_levels);
    }

    // part 1
    bool part1 = false;
    if (part1) {
        std::cout<<"Part 1! "<<risk_levels.size()<<" and "<<risk_levels.at(0).size()<<std::endl;
        astar(risk_levels);

        // Reset the old risk levels/costs.
        for (auto& risk_pt_row : risk_levels) {
            for (auto& risk_pt : risk_pt_row) {
                risk_pt.cost_for_path_to_node = std::numeric_limits<double>::max();
                risk_pt.parent = nullptr;
                risk_pt.f_score = std::numeric_limits<double>::max();
            }
        }
    }

    const std::size_t starting_x = risk_levels.size();
    const std::size_t starting_y = risk_levels.at(0).size();

    // grow risk level map.
    // First grow the rows out horizontally.
    for (int round = 1; round < 5; round++) {
        for (std::size_t i = 0; i < starting_x; i++) {
            // for each row, grow it horizontally first.
            std::vector<RiskPoint> next_addition;
            for (std::size_t j = 0; j < starting_y; j++) {
                int last_y = (round - 1) * starting_y + j;
                double new_risk = risk_levels.at(i).at(last_y).risk + 1;
                if (new_risk > 9) new_risk = 1;
                int new_x = i;
                int new_y = round * starting_y + j;
                next_addition.push_back(RiskPoint(new_x, new_y, new_risk));
            }
            risk_levels.at(i).insert( risk_levels.at(i).end(), next_addition.begin(), next_addition.end() );
        }
    }
    // Next grow the rows down vertically.
    for (int round = 1; round < 5; round++) {
        for (std::size_t i = 0; i < starting_x; i++) {
            std::vector<RiskPoint> next_additions;
            for (std::size_t j = 0; j < risk_levels.at(i).size(); j++) {
                int last_x = (round - 1) * starting_x + i;
                double new_risk = risk_levels.at(last_x).at(j).risk + 1;
                if (new_risk > 9.1) new_risk = 1;
                int new_x = round * starting_x + i;
                int new_y = j;
                next_additions.push_back(RiskPoint(new_x, new_y, new_risk));
            }
            risk_levels.push_back(next_additions);
        }
    }

    if (false) {
        for (std::size_t i = 0; i < risk_levels.size(); i++) {
            std::string x = "";
            for (std::size_t j = 0; j < risk_levels.at(0).size(); j++) {
                x+= std::to_string((int)risk_levels.at(i).at(j).risk);
                // x+= "  ";
            }
            std::cout<<x<<std::endl;
        }
    }

    std::cout<<"part 2 sizes: "<<risk_levels.size()<<" and "<<risk_levels.at(0).size()<<std::endl;
    astar(risk_levels);
}

double apply_type_to_nums(std::vector<double>& nums, int type) {
    if (type == 0) {
        // sum.
        double s = 0;
        for (const double& i : nums) s+= i;
        return s;
    } else if (type == 1) {
        // product
        double s = 1;
        for (const double& i : nums) s *= i;
        return s;
    } else if (type == 2) {
        // min
        return *std::min_element(nums.begin(), nums.end());
    } else if (type == 3) {
        // max.
        return *std::max_element(nums.begin(), nums.end());
    } else if (type == 5) {
        return (nums.at(0) > nums.at(1)) ? 1 : 0;
    } else if (type == 6) {
        return (nums.at(0) < nums.at(1)) ? 1 : 0;
    } else if (type == 7) {
        return (nums.at(0) == nums.at(1)) ? 1 : 0;
    }
    return 0;
}

double handle_subpackets(std::string& full_binary, int& version_sum) {
    std::string version = full_binary.substr(0, 3);
    int version_i = std::stoi(version, 0, 2);
    version_sum += version_i;
    std::string type_id = full_binary.substr(3, 3);
    int type_id_i = std::stoi(type_id, 0, 2);
    full_binary.erase(full_binary.begin(), full_binary.begin()+6);
    if (type_id_i == 4) {
        std::string res = "";
        bool keep_looking = true;
        while (keep_looking) {
            if (full_binary[0] == '0') keep_looking = false;
            res += full_binary.substr(1, 4);
            full_binary.erase(full_binary.begin(), full_binary.begin()+5);
        }
        return std::stoll(res, 0, 2);
    } else {
        // operator packet, so check the next value for size.
        char first_bin = full_binary[0];
        full_binary.erase(full_binary.begin(), full_binary.begin()+1);
        std::vector<double> results;
        if (first_bin == '0') {
            // next 15 bits represents the lenght of total subpackets.
            std::string lengthh = full_binary.substr(0, 15);
            long d_length = std::stoll(lengthh, 0, 2);
            full_binary.erase(full_binary.begin(), full_binary.begin()+15);
            // the full binary string is now full of multiple packets that make up the length of the string.
            // get the result by recursing on this new substring of subpackets. The handle_subpackets function
            // mutates the full_binary string and will continue to truncate it until it is complete, so loop
            // over this recursing on each subpacket.
            std::size_t curr_size = full_binary.size();
            while (full_binary.size() > (curr_size - d_length)) {
                results.push_back(handle_subpackets(full_binary, version_sum));
            }
        } else {
            // next 11 bits represent the number of subpackets.
            std::string num_packets_s = full_binary.substr(0,11);
            double num_packets = std::stoll(num_packets_s, 0, 2);
            full_binary.erase(full_binary.begin(), full_binary.begin()+11);
            // There are num_packets of packets, so iterate that many times with recursion. The handle_subpackets
            // function mutates the full_binary string and will continue to truncate it until complete.
            for (int i = 0; i < num_packets; i++) {
                results.push_back(handle_subpackets(full_binary, version_sum));
            }
        }
        // Now apply the operator from the current packet to each of the subpacket results.
        double op_res = apply_type_to_nums(results, type_id_i);
        return op_res;
    }
    return 0;
}

void day16() {
    std::ifstream input( "./day16data.txt" );
    if(!input.is_open()){
    std::cerr <<  "file cannot be opened";
    }
    if (!input){
    std::cerr << "errors in file";
    }
    std::string line = "";
    std::map<char, std::string> hex_to_bin = {{'0',"0000"}, {'1', "0001"}, {'2',"0010"}, {'3', "0011"}, {'4',"0100"}, {'5',"0101"}, {'6',"0110"}, {'7',"0111"},
                                    {'8',"1000"}, {'9',"1001"},{'A',"1010"}, {'B',"1011"},{'C',"1100"},{'D',"1101"},{'E',"1110"},{'F',"1111"}};
    std::string full_binary = "";
    while (std::getline(input, line)) {
        for (const char& c : line) {
            full_binary+=hex_to_bin[c];
        }
    }
    int version_sum = 0;
    auto res = handle_subpackets(full_binary, version_sum);
    std::cout<<"Version sum: "<<version_sum<<std::endl;
    std::cout<<"Operations result: "<<std::fixed<<res<<std::endl;
}

bool is_in_target(const int& x_min, const int& x_max, const int& y_min, const int& y_max, const int& x, const int& y) {
    // y is negative, so signs are different for the check.
    return ((x >= x_min && x <= x_max) && (y <= y_max && y >= y_min));
}

bool has_overshot_target(const int& x_max, const int& y_min, const int& y_max, const int& x, const int& y) {
    bool y_check = (y_max < y_min) ? y < y_max : y < y_min;
    return x > x_max || y_check;
}

int simulate_until_in_or_past_target(const int& x_min, const int& x_max, const int& y_min, const int& y_max,
                                      const int& x_change, const int& y_change, const int& x_start, const int& y_start) {
    // x increases by x_spec (xspec = xchange+1 if xchange < 0 or xchange-1 if xchange > 0 or xchange if xchange==0)
    // y increases by ychange - 1
    int current_x = x_start;
    int current_y = y_start;
    int curr_x_vel = x_change;
    int curr_y_vel = y_change;
    int highest_y_seen = 0;
    bool keep_going = true;
    // std::cout<<"STARTING WITH "<<curr_x_vel<<" "<<curr_y_vel<<std::endl;
    while (keep_going) {
        // Update the point.
        current_x += curr_x_vel;
        current_y += curr_y_vel;
        // std::cout<<"Current "<<current_x<<" "<<current_y<<std::endl;

        // best y value.
        if (current_y > highest_y_seen) highest_y_seen = current_y;

        // check if we are in the goal!
        if (is_in_target(x_min, x_max, y_min, y_max, current_x, current_y)) {
            // woo! return highest y
            return highest_y_seen;
        }

        // check if we have overshot.
        if (has_overshot_target(x_max, y_min, y_max, current_x, current_y)) {
            return -1000;
        }

        // Update the velocities.
        if (curr_x_vel > 0) curr_x_vel -= 1;
        else if (curr_x_vel < 0) curr_x_vel += 1;
        curr_y_vel -= 1;
    }
    return -1000;
}

void day17() {
    // puzzle data: target area: x=81..129, y=-150..-108
    int x_max = 129;
    int x_min = 81;
    int y_min = -150;
    int y_max = -108;
    // int x_max = 30;
    // int x_min = 20;
    // int y_min = -10;
    // int y_max = -5;
    int x_start = 0;
    int y_start = 0;


    // can bound x/y values to the top of the range because otherwise it will instantly overshoot in the
    // very first step.
    int highest_y = -10000;
    int count_valid_vels = 0;
    for (int xv = 0; xv < x_max + 1; xv++) {
        for (int yv = 0; yv < std::abs(y_min) + 1; yv++) {
            int y_val = simulate_until_in_or_past_target(x_min, x_max, y_min, y_max, xv, yv, x_start, y_start);
            if (y_val > highest_y) highest_y = y_val;
            if (y_val != -1000) count_valid_vels++;
        }
        for (int yv = -1; yv > y_min - 1; yv--) {
            int y_val = simulate_until_in_or_past_target(x_min, x_max, y_min, y_max, xv, yv, x_start, y_start);
            if (y_val > highest_y) highest_y = y_val;
            if (y_val != -1000) count_valid_vels++;
        }

    }
    std::cout<<"Highest y val "<<highest_y<<std::endl;
    std::cout<<"Count of valid velocity pairs "<<count_valid_vels<<std::endl;
}

struct Point3d {
    double x;
    double y;
    double z;

    Point3d(double x_in, double y_in, double z_in) {
        x = x_in;
        y = y_in;
        z = z_in;
    }

};


struct Quat {
    double w;
    double x;
    double y;
    double z;

    Quat(double win, double xin, double yin, double zin) {
        w=win;
        x=xin;
        y=yin;
        z=zin;
    }

    Quat operator*(const Quat& other) const {
        return Quat(
            w * other.w - x * other.x - y * other.y - z * other.z,
            w * other.x + x * other.w + y * other.z - z * other.y,
            w * other.y - x * other.z + y * other.w + z * other.x,
            w * other.z + x * other.y - y * other.x + z * other.w
        );
    }

    Point3d rotate_point(const Point3d& pt) const {
        Quat q = Quat(0, pt.x, pt.y, pt.z);
        Quat inv = Quat(w, -x, -y, -z);
        Quat t1 = q * inv;
        Quat t2 = *this * t1;
        return Point3d(t2.x, t2.y, t2.z);
    }
};

Quat normalize(const Quat& q) {
    double len = std::sqrt((q.w*q.w) + q.x*q.x + q.y*q.y + q.z*q.z);
    if (len < 0.0000001) {
        return Quat(1, 0, 0, 0);
    }
    return Quat(q.w/len, q.x/len, q.y/len, q.z/len);
}

struct PointCloud {
    std::vector<Point3d> points;

    PointCloud apply_rotation(const Quat& rot) const {
        PointCloud rotated_cloud;
        for (const Point3d& pt : points) {
            rotated_cloud.points.push_back(rot.rotate_point(pt));
        }
        return rotated_cloud;
    }

    std::vector<std::pair<int, double>> n_nearest_neighbor_to(const Point3d& input, const int& n) const {
        std::set<int> used_indices;
        std::vector<std::pair<int, double>> res;
        for (int i = 0; i < n; i++) {
            double l1dist = std::numeric_limits<double>::max();
            int best_ind = 0;
            for (std::size_t i = 0; i < points.size(); i++) {
                if (used_indices.count(i) == 1) continue;
                const Point3d& point = points.at(i);
                double this_dist = std::abs(input.x - point.x) + std::abs(input.y - point.y) + std::abs(input.z - point.z);
                if (this_dist < l1dist) {
                    best_ind = i;
                    l1dist = this_dist;
                }
            }
            res.push_back(std::make_pair(best_ind, l1dist));
            used_indices.insert(best_ind);
        }
        return res;
    }
};


// if it is in the same axis system, then 12 beacons will overlap. meaning that there will be at least 12 points with
// the same distance of separation.
// could try rotating then checking how many are nearby
// looking for pairs of scanner which share 12 beacons.
// define everything relative to scanner 0.


void day19() {
    std::ifstream input( "./day19data.txt" );
    if(!input.is_open()){
    std::cerr <<  "file cannot be opened";
    }
    if (!input){
    std::cerr << "errors in file";
    }
    std::string line = "";

    std::vector<PointCloud> scanner_points;
    PointCloud current;
    std::size_t scanner_count = 0;
    while (std::getline(input, line)) {
        size_t n = std::count(line.begin(), line.end(), '-');
        if (n > 3){
            scanner_count+=1;
            // new scanner.
            current = PointCloud();
        } else if (line.empty()) {
            scanner_points.push_back(current);
        }else {
            // add a point to the current scanner.
            auto pt_vec = split_string_to_int(line, ',');
            current.points.push_back(Point3d(pt_vec.at(0), pt_vec.at(1), pt_vec.at(2)));
        }
    }
    if (scanner_count > scanner_points.size()) {
        scanner_points.push_back(current);
    }
    std::cout<<"Scanners "<<scanner_count<<" vs "<<scanner_points.size()<<std::endl;

    std::vector<Quat> rotations = {Quat(1,0,0,0), Quat(1, 0, 1, 0), Quat(0,0,1,0), Quat(1,0,-1,0),
                                   Quat(1, 0, 0, 1), Quat(1, 1, 1, 1), Quat(0, 1, 1, 0), Quat(1, -1, -1, 1),
                                   Quat(1, 0, 0, -1), Quat(1, -1, 1, -1), Quat(0, -1, 1, 0), Quat(1, 1, -1, -1),
                                   Quat(1, 1, 0, 0), Quat(1, 1, 1, -1), Quat(0, 0, 1, -1), Quat(1, 1, -1, 1),
                                   Quat(0, 1, 0, 0), Quat(0, 1, 0, -1), Quat(0, 0, 0, 1), Quat(0, 1, 0, 1),
                                   Quat(1, -1, 0, 0), Quat(0, -1, 1, 1), Quat(0, 0, 1, 1), Quat(1, -1, -1, -1)};
    std::vector<Quat> norm_rotations;
    for (Quat& q : rotations) {
        norm_rotations.push_back(normalize(q));
    }

    // look at pairs. First try to find the first pair.
    const PointCloud& main_cloud = scanner_points.at(0);
    double found_match_dist = -1;
    std::size_t large_num = 12;
    for (std::size_t i = 1; i < scanner_points.size(); i++) {
        // For each rotation, compute nearest negihbors.
        const PointCloud& cloud = scanner_points.at(i);
        for (Quat& rot : norm_rotations) {
            std::cout<<"ROTATION BEING APPLIED "<<rot.w<<" "<<rot.x<<" "<<rot.y<<" "<<rot.z<<std::endl;
            PointCloud rotated_cloud = cloud.apply_rotation(rot);
            // Get nearest neighbors for each rotated point to the main cloud.
            // The pair is (index in main cloud, index in current cloud).
            std::map<double, std::vector<std::pair<int, int>>> distance_and_pts;
            std::map<int, std::set<int>> main_to_nearby_in_cloud;
            for (std::size_t k = 0; k < rotated_cloud.points.size(); k++) {
                const Point3d& rot_pt = rotated_cloud.points.at(k);
                if (rot_pt.x + 68 == 390) {
                std::cout<<"Transformed rot pot: "<<rot_pt.x + 68<<" "<<rot_pt.y + 1246<< " "<<rot_pt.z + 43<<std::endl;

                }
                // TODO: nearest neighbor may be wrong here.
                // try all possible shifts for a specific rotation and if there are > 12 points with distnace of
                // 0 to the nearest neighbor in main, then this is the tf!
                std::vector<std::pair<int, double>> res = main_cloud.n_nearest_neighbor_to(rot_pt, 12);
                for (const auto& n : res) {
                    if (main_to_nearby_in_cloud.count(n.first) > 0) {
                        main_to_nearby_in_cloud.at(n.first).insert((int)k);
                    } else {
                        main_to_nearby_in_cloud[n.first] = {(int)k};
                    }
                }
            }
            // check if there are 12 points which share a large number of nearby points
            int counting_good_cloud_pts = 0;
            for (const auto& cloud_pt_pair : main_to_nearby_in_cloud) {
                if (cloud_pt_pair.second.size() > large_num) counting_good_cloud_pts++;
            }
            std::cout<<"Good cloud pts "<<counting_good_cloud_pts<<std::endl;
            // check if there is a distance with at least 12 points.
            for (const auto& dists_pair : distance_and_pts) {
                // std::cout<<"Dis pairs "<<dists_pair.second.size()<<" with dist "<<dists_pair.first<<std::endl;
                if (dists_pair.second.size() > 12) {
                    found_match_dist = dists_pair.first;
                    break;
                }
            }
            if (found_match_dist != -1) {
                // generate the tf.
                // main_T_cloud.rotation = rot;
                // to get translation, iterate over the points and find the difference from
                // main - cloud.
                for (const auto& pt_pair : distance_and_pts[found_match_dist]) {
                    const Point3d& main_pt = main_cloud.points.at(pt_pair.first);
                    const Point3d& curr_pt = rotated_cloud.points.at(pt_pair.second);
                    std::cout<<"Diff "<<main_pt.x - curr_pt.x<<" " <<main_pt.y - curr_pt.y<<std::endl;
                }
            }
        }
    }

    // accumulate each scanner's beacons realitive to it. Don't include the satelite itself
    // ICP to match each point cloud and determine the tf.

    /**
     * algorithm ICP(M, S)
    θ := θ0
    while not registered:
        X := ∅
        for mi ∊ T(M, θ):
            ŝj := closest point in S to mi
            X := X + ⟨mi, ŝj⟩
        θ := least_squares(X)
    return θ
    Define S = second point set, M = stationary pointset.
    iteratively:
    1. Finds the closest point in S for every point in M
    2. Find the best rigid transformation by solving the least squares problem.
    */

}

int main(int /*argc*/, char** /**argv[]*/)
{
    day17();
}