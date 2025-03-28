#include <iostream>

#include <NTL/ZZ_p.h>
#include <NTL/ZZ.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>

using namespace std;
using namespace NTL;

// Bool_Calculation
namespace Bool_Calculation{
    /////////////////////////////////////////////////////////////////////
    // 基础布尔电路
    bool AND(bool a, bool b) { return a && b; } // AND
    bool O_R(bool a, bool b) { return a || b; } // OR
    bool XOR(bool a, bool b) { return a ^ b; }  // XOR
    bool NOT(bool a) { return !a; }             // NOT

    // 布尔向量初始化：输入一个 Vec<bool> 向量 v1，和 位宽 bitWidth，完成初始化
    void Vec_bool_Initialization(Vec<bool> &v1, int bitWidth){
        v1.SetLength(bitWidth);                                 // 设定长度
        for (int i = 0; i < bitWidth; i++) { v1[i] = 0; }       // 元素置零
    }

    // Half Adder (半加器)：输入两个单比特，输出 sum 和进位 carry
    void HalfAdder(bool a, bool b, bool &sum, bool &carry) {
        sum = XOR(a, b);
        carry = AND(a, b);
    }

    // Half Subtractor（半减器）：输入两个单比特，输出 sum 和借位 borrow_out
    void HalfSubtractor(bool a, bool b, bool &sub, bool &borrow_out) {
        borrow_out = AND(NOT(a), b);    // 计算借位
        sub = XOR(a, b);                // 计算差值   
    }

    // Full Adder (全加器)：输入两个单比特和进位 carry_in，输出 sum 和进位 carry_out
    void FullAdder(bool a, bool b, bool carry_in, bool &sum, bool &carry_out) {
        bool s1, c1, c2;
        HalfAdder(a, b, s1, c1);
        HalfAdder(s1, carry_in, sum, c2);
        carry_out = O_R(c1, c2);
    }

    // Full Subtractor (全减器)：输入两个单比特和借位 borrow_out，输出 sum 和借位 borrow_out
    void FullSubtractor(bool a, bool b, bool borrow_in, bool &sub, bool& borrow_out) {
        bool s1, b1, b2;
        HalfSubtractor(a, b, s1 ,b1);
        HalfSubtractor(s1, borrow_in, sub, b2);
    
        // 计算最终借位：B_out = B1 ∨ B2
        borrow_out = O_R(b1, b2);
    }

    // bool 向量 加法 逻辑：输入两个比特串 v1 v2 和 位宽 bitWidth，输出加法结果 v3 和进位 borrow_out
    void Bool_Add(Vec<bool> v1, Vec<bool> v2, int bitWidth, Vec<bool> &v3, bool& carry_out){
        Vec<bool> temp_sum;             // 定义中间变量 加法结果
        bool temp_carry = 0;            // 定义中间变量 储存进位

        Vec_bool_Initialization(temp_sum, bitWidth);        // 初始化，加法长度不超过位宽

        for (int i = 0; i < bitWidth; i++){    // // 循环加法，从低位开始，[0] 是最低位
            FullAdder(v1[i], v2[i], temp_carry, temp_sum[i], temp_carry);    // 加法
        }

        v3 = temp_sum;
        carry_out = temp_carry;

        //cout << "----ADD----" << endl;
        //cout << v1 << endl;
        //cout << v2 << endl;
        //cout << v3 << endl;
        //cout << temp_carry << endl;
        //cout << "----ADD----" << endl;
    }

    // bool 向量 减法 逻辑：输入两个比特串 v1 v2 和 位宽 bitWidth，输出减法结果 v3 和借位 borrow_out
    void Bool_Sub(Vec<bool> v1, Vec<bool> v2, int bitWidth, Vec<bool> &v3, bool& borrow_out){ 
        Vec<bool> temp_sub;             // 定义中间变量 减法结果
        bool temp_borrow = 0;           // 定义中间变量 储存借位
        
        Vec_bool_Initialization(temp_sub, bitWidth);        // 初始化，减法长度不超过位宽

        for (int i = 0; i < bitWidth; i++){    // 循环减法，从低位开始，[0] 是最低位
            FullSubtractor(v1[i], v2[i], temp_borrow, temp_sub[i], temp_borrow);    // 减法
        }

        v3 = temp_sub;
        borrow_out = temp_borrow;

        //cout << "----SUB----" << endl;
        //cout << v1 << endl;
        //cout << v2 << endl;
        //cout << v3 << endl;
        //cout << borrow_out << endl;
        //cout << "----SUB----" << endl;
    }

    // bool 比较器：如果 v1 >= v2 返回 1 （使用减法Bool_Sub实现）
    bool BoolComparator(Vec<bool> v1, Vec<bool> v2){
        Vec<bool> v3;
        bool borrow_out;
        int bitWidth = v1.length();

        Bool_Sub(v1, v2 , bitWidth, v3, borrow_out);
        return NOT(borrow_out);     // 如果有余数就是小于

    }

    // bool 向量 除法 逻辑：输入两个比特串 v1 v2 和 位宽 bitWidth，输出除法结果 v3 和余数 remainder
    void Bool_Div(Vec<bool> v1, Vec<bool> v2, int bitWidth, Vec<bool> &v3, Vec<bool>& remainder){ 
        Vec<bool> temp_Div, temp_remainder;
        
        Vec_bool_Initialization(temp_Div, bitWidth);        // 初始化
        Vec_bool_Initialization(temp_remainder, bitWidth);  // 初始化

        // [0] 是最低位，除法从最高位开始
        for (int i = bitWidth - 1; i >= 0; --i) {
            // 左移 remainder，空出最低位（但是在这里是右移，因为 [0] 是最低位）
            for (int j = bitWidth - 1; j > 0; --j) {
                temp_remainder[j] = temp_remainder[j - 1];
            }
            // 让最低位 = 被除数的最高位
            temp_remainder[0] = v1[i];
    
            // 试探减法 remainder - v2
            Vec<bool> temp_sub;             // 定义中间参数，保存减法结果
            bool temp_borrow;               // 定义中间参数，保存减法借位

            temp_sub.SetLength(bitWidth);   // 定义中间参数长度
            Bool_Sub(temp_remainder, v2, bitWidth, temp_sub, temp_borrow); // 计算减法 remainder - v2
    
            // 如果 remainder >= v2，即没有借位，则更新
            temp_Div[i] = AND(NOT(temp_borrow), 1); 
            for (int k = 0; k < bitWidth; ++k) {
                temp_remainder[k] = XOR(AND(NOT(temp_borrow), temp_sub[k]), AND(temp_borrow, temp_remainder[k]));
            } 
        }
        
        v3 = temp_Div;
        remainder = temp_remainder;
    
        //cout << "----DIV----" << endl;
        //cout << v1 << endl;
        //cout << v2 << endl;
        //cout << v3 << endl;
        //cout << remainder << endl;
        //cout << "----DIV----" << endl;
    }

    // bool 向量 乘法 逻辑：输入两个比特串 v1 v2 和 位宽 bitWidth，输出减法结果 v3(长 2 * bitWidth)
    void BoolMultiplier(Vec<bool>& v1, Vec<bool>& v2, int bitWidth, Vec<bool>& v3) {
        Vec<bool> temp_mul, shiftedMultiplicand;    // 定义 中间变量 除法结果 左移量
        bool carry_out;                             // 定义 进位
        int result_Width = 2 * bitWidth;            // 结果最大长度为 2 bitWidth

        Vec_bool_Initialization(v3, result_Width);                  // 初始化，结果长度为 2 bitWidth
        Vec_bool_Initialization(temp_mul, result_Width);            // 初始化，结果长度为 2 bitWidth
        Vec_bool_Initialization(shiftedMultiplicand, result_Width); // 初始化，结果长度为 2 bitWidth

        for (int i = 0; i < bitWidth; i++) {
            Vec_bool_Initialization(shiftedMultiplicand, result_Width);     // 初始化，结果长度为 2 bitWidth
            for (int j = 0; j < bitWidth; j++) {
                shiftedMultiplicand[j + i] = AND(v2[i], v1[j]);             // 左移 i 位，并加入掩码，只有当 v2[i] == 1 时才进行累加，代替 if 语句
            }
    
            temp_mul = v3; // 备份当前 v3
            Bool_Add(temp_mul, shiftedMultiplicand, result_Width, v3, carry_out); // 加法
            
        }

        //cout << "----MUL----" << endl;
        //cout << v1 << endl;
        //cout << v2 << endl;
        //cout << v3 << endl;
        //cout << "----MUL----" << endl;
    }

    // bool 向量 mod 逻辑：输入两个比特串 v1 v2 和 结果的位宽 bitWidth，输出 mod 数结果(bitWidth，输出)
    void Bool_Mod(Vec<bool>& v1, Vec<bool>& v2, int bitWidth, Vec<bool>& v_mod){
        Vec<bool> v3, temp_v2, temp_mod;
        Vec_bool_Initialization(temp_v2, v1.length()); // 初始化，结果长度为 v1.length()

        // 让 v2 扩展到 v1 的长度，防止除法越界
        for (int i = 0; i < bitWidth; i++){
            temp_v2[i] = v2[i];
        }
        for (int i = bitWidth; i < v1.length(); i++){
            temp_v2[i] = 0;
        }

        // 进行除法，余数 就是答案
        Bool_Div(v1, temp_v2, v1.length(), v3, temp_mod);   

        // 舍弃需要位宽之后的位
        for (int i = 0; i < bitWidth; i++){
            v_mod[i] = temp_mod[i];
        }
        
        //cout << "----MOD----" << endl;
        //cout << v1 << endl;
        //cout << v2 << endl;
        //cout << v_mod << endl;
        //cout << "----MOD----" << endl;
    }

    // ZZ_p 的布尔矩阵乘法：M3 = M1 * M2
    Mat<Vec<bool>> bool_Matrix_Mul(Mat<Vec<bool>> M1, Mat<Vec<bool>> M2, Vec<bool> q){
        Mat<Vec<bool>> M3;                  // 结果矩阵
        M3.SetDims(M1.NumRows(), M2.NumCols());

        int bitWidth = M1[0][0].length();   // 位宽，也就是向量长度
        
        // M3 归零
        for (int i = 0; i < M3.NumRows(); i++) {
            for (int j = 0; j < M3.NumCols(); j++) {
                M3[i][j].SetLength(bitWidth);
                for (int k = 0; k < bitWidth; k++){
                    M3[i][j][k] = 0;
                }
            }
        }
        // 乘法
        for (int i = 0; i < M1.NumRows(); i++) {
            for (int j = 0; j < M2.NumCols(); j++) {
                for (int k = 0; k < M1.NumCols(); k++) {
                    Vec<bool> temp, temp_mod, temp_sum;
                    bool temp_carry;

                    temp_mod.SetLength(bitWidth);

                    BoolMultiplier(M1[i][k], M2[k][j], bitWidth, temp);         // 乘法
                    Bool_Mod(temp, q, bitWidth, temp_mod);                          // mod
                    Bool_Add(M3[i][j], temp_mod, bitWidth, temp_sum, temp_carry);    // 加法
                    
                    Vec<bool> temp_sum_carry;                   
                    temp_sum_carry.SetLength(bitWidth + 1); 

                    for (int i = 0; i < bitWidth; i++){
                        temp_sum_carry[i] = temp_sum[i];
                    }
                    temp_sum_carry[bitWidth] = temp_carry;          // 加上进位

                    Bool_Mod(temp_sum_carry, q, bitWidth, M3[i][j]);    // 再 mod 一次
                }
            }
        }

        return M3;
    }
}
using namespace Bool_Calculation;

int main() {
    cout << "#------------- Start -------------#\n" << endl;
    
    cout << "#----------------------- End of Program -----------------------#" << endl;

    return 0;
}
