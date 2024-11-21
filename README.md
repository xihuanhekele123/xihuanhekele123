#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <algorithm> // 引入算法库以使用min函数

using namespace std;

// 1初始化所有参数
const int ndivx = 300; // x方向质点总数
const int ndivy = 88; // y方向质点总数
const int nbnd = 3;
const int totnode = (ndivx + 2 * nbnd) * (ndivy + nbnd);
const int tot = ndivx * ndivy; // 内部总质点数

const int nt = 1000; // 总时间步数
const int maxfam = 100;
const double length = 6e-1; // 板的长度
const double width = 176e-3; // 板的宽度
const double dx = length / ndivx; // 质点间距
const double dy = width / ndivy;
const double delta = 3.015 * dx; // 领域半径
const double dens = 7850.0; // 密度
const double emod = 21e10; // 弹性模量（杨氏）
const double pratio = 1.0 / 3.0; // 泊松比
const double alpha = 23e-6; // 热膨胀系数
const double dtemp = 0.0; // 温度变化
const double area = dx * dx; // 横截面积
const double vol = area * dx; // 单个质点体积
const double bc = 9.0 * emod / (M_PI * dx * (delta * delta * delta)); // 粘结常数
const double sedload1 = 9.0 / 16.0 * emod * 1e-6; // 第一加载条件下的应变能密度
const double sedload2 = 9.0 / 16.0 * emod * 1e-6; // 第二加载条件下的应变能密度
const double dt = 1.0; // 时间步长
const double totime = nt * dt; // 总时间
const double current_time = 0.0; // 当前时间
const double a = 65e-4;
const double b = 45e-4;
const double P = 69800;
const double appres = 200e6; // 拉伸荷载
const double radij = dx / 2.0; // 材料点半径
const double fc = 0.3;
const double scr0 = 0.005;
const double A1 = 426.0, m1 = 2.77;
double lambda = 1.0; // 初始寿命
double currenttime = 0.0;
double dforce1;
double dforce2;
double Ni_min=1e10;
int totint = 0;
int totintleft=0;
int totintright=0;
int totintbottom =0;

// 动态数组存储质点的信息
vector<vector<double>> coord(totnode, vector<double>(2, 0.0)); // 质点坐标
vector<vector<double>> pforce(totnode, vector<double>(2, 0.0)); // 质点受力
vector<vector<double>> pforceold(totnode, vector<double>(2, 0.0)); // 上一时间步的质点受力
vector<vector<double>> bforce(totnode, vector<double>(2, 0.0)); // 体力
vector<vector<double>> stendens(totnode, vector<double>(2, 0.0)); // 应变能
vector<vector<double>> fncst(totnode, vector<double>(2, 1.0)); // 表面修正因子
vector<vector<double>> disp(totnode, vector<double>(2, 0.0)); // 位移
vector<vector<double>> andisp(totnode, vector<double>(2, 0.0));
vector<vector<double>> vel(totnode, vector<double>(2, 0.0)); // 速度
vector<vector<double>> velhalfold(totnode, vector<double>(2, 0.0)); // 上一时间步的半步速度
vector<vector<double>> velhalf(totnode, vector<double>(2, 0.0)); // 当前半步速度
vector<vector<double>> acc(totnode, vector<double>(2, 0.0)); // 加速度
vector<vector<double>> massvec(totnode, vector<double>(2, 0.0)); // 质量向量
//vector<vector<double>> smax(totnode, vector<double>(2, 0.0)); //
//vector<vector<double>> smin(totnode, vector<double>(2, 0.0)); //
std::vector<double> enddisp(nt, 0.0);
std::vector<double> endtime(nt, 0.0);
std::vector<double> dmg(totnode, 0.0);
std::vector<std::vector<int>> fail(totnode, std::vector<int>(maxfam, 1));
std::vector<std::vector<double>> remainingLife(totnode, std::vector<double>(maxfam, 1));

double max_dmg=0;

//vector<double> Ni_min(totnode, 0.0);
//vector<vector<double>> remainingLife(totnode, vector<double>(2, 1.0)); // 每个键的剩余寿命
//2函数
// 计算两点间距离的函数
double distance(const vector<double>& a, const vector<double>& b) {
    return sqrt(pow(a[0] - b[0], 2) + pow(a[1] - b[1], 2));
}

// 计算角度的函数
double calculateTheta(double dy, double dx) {
    if (abs(dx) < 1e-10) return M_PI / 2;
    if (abs(dy) < 1e-10) return 0;
    return atan2(dy, dx);
}

// 计算键的剩余寿命
double calculate_remaining_life(int i, int j, double smax, double smin, const double A1, const double m1, double Ni,double rL) {
    if (rL > 0) {
        return rL - Ni * A1 * pow(abs(smax-smin), m1);
    } else {
        return 0;
    }
}
// 计算Ni
double calculate_Ni(int i, int j, double smax, double smin, const double A1, const double m1,double rL) {
    if (rL > 0) {
        return (rL) / (A1 * pow(abs(smax-smin), m1));
    }else{
        return 0;
    }
}

// 计算法向接触应力
double normal_stress(double x, double y) {
    if (x >= -a && x <= a && y >= 0.5 * width - dx && y <= 0.5 * width) {
        return 3 * P / (2 * M_PI * a * b) * sqrt(1 - (x * x) / (a * a));
    } else {
        return 0;
    }
}

int main() {
    int nnum = 0;
    // 3设置质点坐标
    for (int i = 0; i < ndivx; ++i) {
        for (int j = 0; j < ndivy; ++j) {
            double coordx = -0.5 * length + (dx / 2.0) + (i) * dx;
            double coordy = -0.5 * width + (dx / 2.0) + (j) * dx;
            coord[nnum][0] = coordx;
            coord[nnum][1] = coordy;
            nnum++;
        }
    }
    totint = nnum;
    // 左边界
    for (int i = 0; i < nbnd; ++i) {
        for (int j = 0; j < ndivy; ++j) {
            double coordx = -0.5 * length - 3 * dx + (dx / 2.0) + (i) * dx;
            double coordy = -0.5 * width + (dx / 2.0) + (j) * dx;
            coord[nnum][0] = coordx;
            coord[nnum][1] = coordy;
            nnum++;
        }
    }
    totintleft=nnum;
    // 右边界
    for (int i = 0; i < nbnd; ++i) {
        for (int j = 0; j < ndivy; ++j) {
            double coordx = 0.5 * length + (dx / 2.0) + (i) * dx;
            double coordy = -0.5 * width + (dx / 2.0) + (j) * dx;
            coord[nnum][0] = coordx;
            coord[nnum][1] = coordy;
            nnum++;
        }
    }
    totintright = nnum;
    // 下边界
    for (int i = 0; i < (ndivx + 2 * nbnd); ++i) {
        for (int j = 0; j < nbnd; ++j) {
            double coordx = -0.5 * length - 3 * dx + (dx / 2.0) + (i) * dx;
            double coordy = -0.5 * width - 3 * dx + (dx / 2.0) + (j) * dx;
            coord[nnum][0] = coordx;
            coord[nnum][1] = coordy;
            nnum++;
        }
    }
    totintbottom=nnum;

    // 4动态数组存储近场范围内质点信息
    std::vector<int> numfam(totnode, 0); // 每个质点在近场范围内的其他点
    std::vector<int> pointfam(totnode, 0); // 进场范围成员索引数组
    std::vector<int> nodefam(10000000, 0); // 近场范围成员数组

    for (int i = 0; i < totnode; ++i){
        if (i == 0) {
            pointfam[i] = 0;
        } else {
            pointfam[i] = pointfam[i - 1] + numfam[i - 1];
        }
        for (int j = 0; j < totnode; ++j) {
            double idist = distance(coord[i], coord[j]);
            if (i != j && idist <= delta) {
                numfam[i]++;
                nodefam[pointfam[i] + numfam[i] - 1] = j;
            }
        }
    }

    // 6计算表面修正因子
    // X方向
    for (int i = 0; i < totnode; ++i) {
        disp[i][0] = 0.001 * coord[i][0];
        disp[i][1] = 0.0;
    }
    for (int i = 0; i < totnode; ++i) {
        stendens[i][0] = 0.0;
        for (int j = 0; j < numfam[i]; ++j) {
            int cnode = nodefam[pointfam[i] + j];
            double idist = distance(coord[i], coord[cnode]);
            double nlength = distance({coord[cnode][0] + disp[cnode][0], coord[cnode][1] + disp[cnode][1]},
                                      {coord[i][0] + disp[i][0], coord[i][1] + disp[i][1]});
            double fac = (idist <= delta - radij) ? 1.0 : ((idist <= delta + radij) ? (delta + radij - idist) /
                                                                                      (2.0 * radij) : 0.0);
            stendens[i][0] += 0.5 * 0.5 * bc * pow((nlength - idist) / idist, 2) * idist * vol * fac;
        }
        fncst[i][0] = sedload1 / stendens[i][0];
    }

// Y方向
    for (int i = 0; i < totnode; ++i) {
        disp[i][0] = 0.0;
        disp[i][1] = 0.001 * coord[i][1];
    }
    for (int i = 0; i < totnode; ++i) {
        stendens[i][1] = 0.0;
        for (int j = 0; j < numfam[i]; ++j) {
            int cnode = nodefam[pointfam[i] + j];
            double idist = distance(coord[i], coord[cnode]);
            double nlength = distance({coord[cnode][0] + disp[cnode][0], coord[cnode][1] + disp[cnode][1]},
                                      {coord[i][0] + disp[i][0], coord[i][1] + disp[i][1]});
            double fac = (idist <= delta - radij) ? 1.0 : ((idist <= delta + radij) ? (delta + radij - idist) /
                                                                                      (2.0 * radij) : 0.0);
            stendens[i][1] += 0.5 * 0.5 * bc * pow((nlength - idist) / idist, 2) * idist * vol * fac;
        }
        fncst[i][1] = sedload2 / stendens[i][1];
    }

// 7初始化速度和位移
    for (int i = 0; i < totnode; ++i) {
        vel[i][0] = 0.0;
        disp[i][0] = 0.0;
        vel[i][1] = 0.0;
        disp[i][1] = 0.0;
    }

// 8计算稳定质量向量
    for (int i = 0; i < totnode; ++i) {
        massvec[i][0] = 0.25 * dt * dt * (M_PI * delta * delta * dx) * bc / dx * 5.0;
        massvec[i][1] = 0.25 * dt * dt * (M_PI * delta * delta * dx) * bc / dx * 5.0;
    }

    // 9应用加载
    for (int i = 0; i < totnode; ++i) {
        const double appres = normal_stress(coord[i][0], coord[i][1]); // 荷载
        bforce[i][1] = -appres / dx;//
        bforce[i][0] =  fc*appres / dx;
        if (bforce[i][1] != 0) {
            std::cout << bforce[i][1] << std::endl;
        }
        if (bforce[i][0] != 0) {
            std::cout << bforce[i][0] << std::endl;
        }
    }

// 文件写入
    ofstream outFile("zzh_coord_disp_damage.txt");
    bool terminate = false;
    while (!terminate) {
        std::vector<std::vector<double>> smax(totnode, std::vector<double>(maxfam, 0));
        std::vector<std::vector<double>> smin(totnode, std::vector<double>(maxfam, 0));
        std::vector<std::vector<double>> Ni(totnode, std::vector<double>(maxfam, 100000000));
        std::vector<std::vector<double>> eps(totnode, std::vector<double>(maxfam,0));
        // 10时间积分
        for (int tt = 0; tt < nt; ++tt) {
            std::cout << "tt = " << tt << std::endl;
            currenttime = tt * dt;
//边界约束
            for (int i = 0; i < totnode; ++i) {
                if (coord[i][1] < -0.5 * width) {
                    vel[i][1] = 0;
                    disp[i][1] = 0;
                }
            }
            for (int i = 0; i < totnode; ++i) {
                if (coord[i][0] < -0.5 * length) {
                    vel[i][0] = 0;
                    disp[i][0] = 0;
                }
            }
            for (int i = 0; i < totnode; ++i) {
                if (coord[i][0] > 0.5 * length) {
                    vel[i][0] = 0;
                    disp[i][0] = 0;
                }
            }

            double max_eps = -1.0;
            int max_eps_i = -1;  // 记录最大 eps 对应的 i
            int max_eps_j = -1;  // 记录最大 eps 对应的 j
            for (int i = 0; i < tot; ++i) {
                double dmgpar1 = 0.0;
                double dmgpar2 = 0.0;
                pforce[i][0] = 0.0;
                pforce[i][1] = 0.0;
                for (int j = 0; j < numfam[i]; ++j) {
                    //smax[i][j] = -numeric_limits<double>::max();
                    //smin[i][j] = numeric_limits<double>::max();
                    int cnode = nodefam[pointfam[i] + j];
                    double idist = distance(coord[i], coord[cnode]);
                    double nlength = distance({coord[cnode][0] + disp[cnode][0], coord[cnode][1] + disp[cnode][1]},
                                              {coord[i][0] + disp[i][0], coord[i][1] + disp[i][1]});
                    double strain = abs((nlength - idist) / idist);
                    smax[i][j] = max(smax[i][j], strain);
                    double s_max = smax[i][j];
                    if (smin[i][j] == 0) {
                        smin[i][j] = strain;
                    } else if (strain != 0 && strain < smin[i][j]) {
                        smin[i][j] = strain;
                    }
                    double s_min = smin[i][j];
                    double fac = (idist <= delta - dx / 2.0) ? 1.0 : (idist <= delta + dx / 2.0) ?
                                                                     (delta + dx / 2.0 - idist) / (2.0 * dx / 2.0)
                                                                                                 : 0.0;
                    double theta = calculateTheta(coord[cnode][1] - coord[i][1], coord[cnode][0] - coord[i][0]);
                    double scx = (fncst[i][0] + fncst[cnode][0]) / 2.0;
                    double scy = (fncst[i][1] + fncst[cnode][1]) / 2.0;
                    double scr =
                            1.0 / ((cos(theta) * cos(theta) / (scx * scx) + sin(theta) * sin(theta) / (scy * scy)));
                    scr = sqrt(scr);


                    if (fail[i][j] == 1) {
                        dforce1 = bc * ((nlength - idist) / idist - alpha * dtemp) * vol * scr * fac *
                                  (coord[cnode][0] + disp[cnode][0] - coord[i][0] - disp[i][0]) / nlength;
                        dforce2 = bc * ((nlength - idist) / idist - alpha * dtemp) * vol * scr * fac *
                                  (coord[cnode][1] + disp[cnode][1] - coord[i][1] - disp[i][1]) / nlength;
                    } else {
                        dforce1 = 0.0;
                        dforce2 = 0.0;
                    }
                    pforce[i][0] += dforce1;
                    pforce[i][1] += dforce2;
#if 0
                    if (std::abs(strain) > scr0) {
                        if (std::abs(coord[i][1]) >= 0) {
                            fail[i][j] = 0;
                            //Ni_min = 1;
                        }
                    }
#endif
                    eps[i][j]=smax[i][j]-smin[i][j];
                    if (eps[i][j] > max_eps) {
                        max_eps = eps[i][j];      // 更新全局最大 eps
                        max_eps_i = i;            // 记录最大 eps 对应的 i
                        max_eps_j = j;            // 记录最大 eps 对应的 j
                    }
                }
            }
            if (max_eps_i != -1 && max_eps_j != -1) {
                fail[max_eps_i][max_eps_j] = 0;  // 标记为断裂
            }
            // 11自适应动态松弛
            double cn = 0.0, cn1 = 0.0, cn2 = 0.0;
            for (int i = 0; i < tot; ++i) {
                if (velhalfold[i][0] != 0.0) {
                    cn1 -= disp[i][0] * disp[i][0] *
                           (pforce[i][0] / massvec[i][0] - pforceold[i][0] / massvec[i][0]) /
                           (dt * velhalfold[i][0]);
                }
                if (velhalfold[i][1] != 0.0) {
                    cn1 -= disp[i][1] * disp[i][1] *
                           (pforce[i][1] / massvec[i][1] - pforceold[i][1] / massvec[i][1]) /
                           (dt * velhalfold[i][1]);
                }
                cn2 += disp[i][0] * disp[i][0] + disp[i][1] * disp[i][1];
            }

            if (cn2 != 0.0) {
                if ((cn1 / cn2) > 0.0) {
                    cn = 2.0 * sqrt(cn1 / cn2);
                } else {
                    cn = 0.0;
                }
            } else {
                cn = 0.0;
            }

            if (cn > 2.0) {
                cn = 1.9;
            }

            for (int i = 0; i < tot; ++i) {
                if (tt == 0) {
                    velhalf[i][0] = 1.0 * dt / massvec[i][0] * (pforce[i][0] + bforce[i][0]) / 2.0;
                    velhalf[i][1] = 1.0 * dt / massvec[i][1] * (pforce[i][1] + bforce[i][1]) / 2.0;
                } else {
                    velhalf[i][0] = ((2.0 - cn * dt) * velhalfold[i][0] +
                                     2.0 * dt / massvec[i][0] * (pforce[i][0] + bforce[i][0])) /
                                    (2.0 + cn * dt);
                    velhalf[i][1] = ((2.0 - cn * dt) * velhalfold[i][1] +
                                     2.0 * dt / massvec[i][1] * (pforce[i][1] + bforce[i][1])) /
                                    (2.0 + cn * dt);
                }

                vel[i][0] = 0.5 * (velhalfold[i][0] + velhalf[i][0]);
                vel[i][1] = 0.5 * (velhalfold[i][1] + velhalf[i][1]);
                disp[i][0] += velhalf[i][0] * dt;
                disp[i][1] += velhalf[i][1] * dt;

                velhalfold[i][0] = velhalf[i][0];
                velhalfold[i][1] = velhalf[i][1];
                pforceold[i][0] = pforce[i][0];
                pforceold[i][1] = pforce[i][1];
            }
            endtime[tt] = currenttime;
        }

        for (int i = 0; i < tot; ++i) {
            double dmgpar1 = 0.0;
            double dmgpar2 = 0.0;
            double remaining_life;
            for (int j = 0; j < numfam[i]; ++j) {
                if (fail[i][j] == 1) continue;
                int cnode = nodefam[pointfam[i] + j];
                double s_max = smax[i][j];
                double s_min = smin[i][j];
                double idist = distance(coord[i], coord[cnode]);
                double nlength = distance({coord[cnode][0] + disp[cnode][0], coord[cnode][1] + disp[cnode][1]},
                                          {coord[i][0] + disp[i][0], coord[i][1] + disp[i][1]});
                double fac = (idist <= delta - dx / 2.0) ? 1.0 : (idist <= delta + dx / 2.0) ?
                                                                 (delta + dx / 2.0 - idist) / (2.0 * dx / 2.0)
                                                                                             : 0.0;
                double theta = calculateTheta(coord[cnode][1] - coord[i][1], coord[cnode][0] - coord[i][0]);
                double scx = (fncst[i][0] + fncst[cnode][0]) / 2.0;
                double scy = (fncst[i][1] + fncst[cnode][1]) / 2.0;
                double scr = 1.0 / ((cos(theta) * cos(theta) / (scx * scx) + sin(theta) * sin(theta) / (scy * scy)));
                scr = sqrt(scr);
// 13. 计算所有键的Ni值，并找到每个质点的最小Ni值及其索引
//vector<vector<double>> Ni_min(totnode, vector<double>(1,numeric_limits<double>::max())); // Ni_min 存储每个质点的最小Ni值
                vector<pair<int, int>> min_Ni_indices(totnode, make_pair(-1, -1)); // 用于存储最小Ni对应的键索引
                // 计算剩余寿命
                remaining_life = calculate_remaining_life(i, j, s_max, s_min, A1, m1, 1,
                                                          remainingLife[i][j]);
                remainingLife[i][j] = remaining_life;

                // 计算当前键的Ni值
                double epsilon = 1e-10; // 一个非常小的常数
                if (remainingLife[i][j] > 0) {
                    Ni[i][j] = remainingLife[i][j] / (A1 * pow(max(abs(s_max-s_min), epsilon), m1));
                }
                cout << "smax: " << smax[i][j] << ", remainingLife: " << remainingLife[i][j] << ", Ni: " << Ni[i][j] << endl;

                // 如果当前键的Ni小于该质点的最小Ni，则更新最小Ni值和对应的键索引
                if (1 < Ni[i][j] < Ni_min) {
                    Ni_min = Ni[i][j]; // 更新最小Ni值
                    min_Ni_indices[i] = make_pair(i, j); // 更新最小Ni对应的键索引
                     if(Ni[i][j]<=1){
                        Ni_min=1;
                    }
                }
                cout << "smax: " << smax[i][j] << ", remainingLife: " << remainingLife[i][j] << ", Ni_min: "
                     << Ni_min << endl;
                if (Ni_min != 0) {
                    std::cout << ", Ni_min: " << Ni_min << std::endl;
                }

                #if 0
// 14. 找到最小Ni值对应的键并断裂
                if (min_Ni_indices[i].first != -1 && min_Ni_indices[i].second != -1) {
                    // 获取最小Ni对应的键的节点索引
                    int i = min_Ni_indices[i].first;
                    int j = min_Ni_indices[i].second;
                    int cnode = nodefam[pointfam[i] + j];

                    // 将对应的键标记为断裂
                    fail[i][j] = 0;
                    // 输出断裂信息（可选）
                    // cout << "Breaking bond between nodes " << min_Ni_indices[i].first << " and " << cnode << endl;
                }
                #endif
                // 15使用Ni_min更新所有键的remaining_life
                remaining_life = calculate_remaining_life(i, j, s_max, s_min, A1, m1, Ni_min,
                                                          remainingLife[i][j]);
                remainingLife[i][j] = remaining_life;
                //计算损伤指数
                dmgpar1 += fail[i][j] * vol * fac;
                dmgpar2 += vol * fac;
            }
            dmg[i] = 1.0 - dmgpar1 / dmgpar2;
            cout << "dmgpar1: " << dmgpar1 << ", dmgpar2: " << dmgpar2 << endl;
            max_dmg = max(max_dmg, dmg[i]);
        }
        // 变形恢复处理
        for (int i = 0; i < tot; ++i) {
            for (int j = 0; j < numfam[i]; ++j) {
                if (fail[i][j] == 0) continue; // 只处理断裂的键
                int cnode = nodefam[pointfam[i] + j];
                double idist = distance(coord[i], coord[cnode]);
                double nlength = distance({coord[cnode][0] + disp[cnode][0], coord[cnode][1] + disp[cnode][1]},
                                          {coord[i][0] + disp[i][0], coord[i][1] + disp[i][1]});
                double strain = (nlength - idist) / idist;

                // 计算恢复位移，这里假设一个简单的线性恢复
                double recovery_distance = strain * idist * 0.5; // 假设恢复一半的伸长量
                vector<double> recovery_direction = {
                        (coord[cnode][0] - coord[i][0]) / idist,
                        (coord[cnode][1] - coord[i][1]) / idist
                };

                // 更新质点位移
                disp[i][0] -= recovery_distance * recovery_direction[0];
                disp[i][1] -= recovery_distance * recovery_direction[1];
                disp[cnode][0] += recovery_distance * recovery_direction[0];
                disp[cnode][1] += recovery_distance * recovery_direction[1];
            }
        }

        for (int i = 0; i < tot; ++i) {
            //disp[i][2]= sqrt(disp[i][0]*disp[i][0]+disp[i][1]*disp[i][1]);总位移
            outFile << coord[i][0] << " " << coord[i][1] << " " << dmg[i] << std::endl;
            //" " << disp[i][0] << " " << disp[i][1] <<" " << Ni_min <<
        }

        // 17检查是否达到终止条件
        if (max_dmg >= 0.5) {
            cout << "Damage index reached " << max_dmg << ", terminating calculation." << endl;
            terminate = true;
        } else {
            cout << "Damage index is " << max_dmg << ", continuing calculation." << endl;
        }
    }
// 关闭文件
    outFile.close();
    return 0;
}
