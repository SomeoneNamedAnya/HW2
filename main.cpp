#include "./scary.h"
#include <sstream>
#include <fstream>

using namespace std;
int main(int argc, char* argv[]) {
    
    
    // Парсинг параметров с которыми запускается программа
    // В param хранятся строки нужных типов и название файла
    // [0] = --p-type= ; [1] = --v-type=, [2] = --v-flow-type=
    // имя файла filename
    vector<pair<int, pair<int, int>>> param(3);
    string filename;
    for (int i = 1; i < argc; i++) {
        //cout << i << "\n";
        string temp = argv[i];
        //cout << temp << endl;
        string flag, ans;
        int fl = -1;
        for (int j = 0; j < temp.size(); j++) {
            if (fl != -1) {
                ans.push_back(temp[j]);
                continue;
            }
            flag.push_back(temp[j]);
            if (flag == "--p-type") {
                fl = 0;
                j += 1;
            } else if (flag == "--v-type") {
                fl = 1;
                j += 1;
            } else if (flag == "--v-flow-type") {
                fl = 2;
                j += 1;
            }
        }
        if (fl == -1) {
            filename = flag;
        } else {
            if (ans == "DOUBLE") {
                param[fl] = DOUBLE;
            } else if (ans == "float") {
                param[fl] = FLOAT;
            } else {
                string first, sec, th;
                int n1 = 0, m1 = 0;
                int ff = 0;
                for (int j = 0; j < ans.size(); j++) {
                    if (ans[j] == '(') {
                        ff = 1;
                        continue;
                    } else if (ans[j] == ',') {
                        ff = 2;
                        continue;
                    } else if (ans[j] == ')') {
                        break;
                    }

					if(ff == 0) {
						first.push_back(ans[j]);
					} else if (ff == 1) {
						n1 *= 10;
						n1 += ans[j] - '0';
					} else {
						m1 *= 10;
						m1 += ans[j] - '0';
					}
                    
                }
                if (first == "FIXED") {
                    param[fl] = make_pair(3, make_pair(n1, m1));
                } else if (first == "FAST_FIXED") {
                    param[fl] = make_pair(4, make_pair(n1, m1));
                }

            }
        }
        
    }
    cout << "Selected parameters:\n";
    cout << param[0].first << " " << param[0].second.first << " " << param[0].second.second << "\n";
    cout << param[1].first << " " << param[1].second.first << " " << param[1].second.second << "\n";
    cout << param[2].first << " " << param[2].second.first << " " << param[2].second.second << "\n";

    // инстанцирование шаблонов, функция constexpr
    instantiation();
    // Полученные шаблоны лежат в map<pair<pair<int, int>, pair<int, int>, void *> ref_er 
    // Где под void * лежит class Base = 0 и Simulation наследуется от Base (Fixed не подменяется другим классом)
    
    cout << " Количество созданных шаблонов = " << ref_er.size() << "\n";
    
    // Считыване данных из файла
    int input_n, input_m;
    size_t init_Time;

    ifstream file;
    file.open(filename);
    if (!file.is_open()) {
        cout << "Файл " << filename << "не удалось открыть :(\n";
        exit(-1);
    }
    string tempstring;
    getline(file, tempstring);
    istringstream SS(tempstring);
    SS >> input_n >> input_m >> init_Time;
    
    vector<vector<char>>  te(input_n, vector<char> (input_m));
    for (int i = 0; i < input_n; i++) {
       
        getline(file, tempstring);
        cout <<tempstring << endl;
        for (int j = 0; j < input_m; j++) {
            te[i][j] = tempstring[j];
           
        }
    }

    // Проверка есть ли нужный размер поля в size__ (предопределен ли он)
    bool is_size_here = false;
    for (int i= 0; i < sz_size; i++) {
        if (make_pair(input_n, input_m) == size__[i]) {
            is_size_here = true;
            break;
        }
    }
     
    // поиск нужного шаблона класса Simulation, инициализация и старт
    int ind1 = 0, ind2 = 0, ind3 = 0, ind4 = 0;
    int flag_for_start = 0;
    for (auto elem : ref_er) {
        ind1 = elem.first.first.first;
        ind2 = elem.first.first.second;
        ind3 = elem.first.second.first;
        ind4 = elem.first.second.second;
    
        if (!is_size_here) {
            if (sz_size == 0 && ind4 == 1) {
                if (param[0] == types[ind1] && param[1] == types[ind2] &&  param[2] == types[ind3]) {
                    static_cast<Base *>(elem.second) -> init_m(te, init_Time);
                    static_cast<Base *>(elem.second) -> start();
                    cout << "End of Simulation\n";
                    exit(0);
                }
            } else if (ind4 == sz_size) {
                if (param[0] == types[ind1] && param[1] == types[ind2] &&  param[2] == types[ind3]) {
                    static_cast<Base *>(elem.second) -> init_m(te, init_Time);
                    static_cast<Base *>(elem.second) -> start();
                    cout << "End of Simulation\n";
                    exit(0);
                }
            }
            continue;
        } else {
            if (param[0] == types[ind1] && param[1] == types[ind2] &&  param[2] == types[ind3] 
                && input_n == size__[ind4].first && input_m == size__[ind4].second) {
                    
                    static_cast<Base *>(elem.second) -> init(te, init_Time);
                    static_cast<Base *>(elem.second) -> start();
                    cout << "End of Simulation\n";
                    exit(0);
                }
        }
    } 
    cout << "Template is not found...";
    
}
