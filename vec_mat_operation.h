#ifndef _MYFOOD_VEC_OP_H
#define _MYFOOD_VEC_OP_H

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <assert.h>
#include <math.h>
#include <map>

template<typename T>
void vec_add(std::vector<T> &dst_vec, const std::vector<T> &src_vec)
{
    assert(dst_vec.size()==src_vec.size());

    for(std::size_t i=0; i<dst_vec.size(); i++)
    {
        dst_vec[i]+=src_vec[i];
    }
}

template<typename T>
void vec_sub(std::vector<T> &dst_vec, const std::vector<T> &src_vec)
{
    assert(dst_vec.size()==src_vec.size());

    for(std::size_t i=0; i<dst_vec.size(); i++)
    {
        dst_vec[i]-=src_vec[i];
    }
}

template<typename T>
void vec_dot_mult(std::vector<T> &dst_vec, const T &src)
{
    for(std::size_t i=0; i<dst_vec.size(); i++)
    {
        dst_vec[i]*=src;
    }
}

template<typename T>
void vec_dot_div(std::vector<T> &dst_vec, const T &src)
{
    for(std::size_t i=0; i<dst_vec.size(); i++)
    {
        dst_vec[i]/=src;
    }
}

template<typename T>
T vec_fabs_sum(const std::vector<T> &src_vec)
{
    T s=0;

    for(std::size_t i = 0; i<src_vec.size(); i++)
    {
        s+=fabs(src_vec[i]);
    }
    return s;
}

template<typename T>
T vec_dot_product(const std::vector<T> &vec1, const std::vector<T> &vec2)
{
    T p=0;

    assert(vec1.size()==vec2.size());

    for(std::size_t i=0; i<vec1.size(); i++)
    {
        p+=vec1[i]*vec2[i];
    }
    return p;
}

template<typename T>
void output_vec(const std::vector<T> &vec)
{
    for(std::size_t i=0; i<vec.size(); i++)
    {
        std::cout << vec[i] << " ";
    }
    std::cout << std::endl;
}

template<typename T>
void vec_mult_to_mat(const std::vector<T> &col_vec, const std::vector<T> &row_vec, std::vector<std::vector<T> > &mat)
{
    mat.clear();

    assert(col_vec.size()==row_vec.size());

    std::vector<T> row(col_vec.size(), 0);

    mat.resize(col_vec.size(), row);

    for(std::size_t r=0; r<col_vec.size(); r++)
    {
        for(std::size_t c=0; c<row_vec.size(); c++)
        {
            mat[r][c]=col_vec[r]*row_vec[c];
        }
    }
}

template<typename T>
void mat_add(std::vector<std::vector<T> > &dst_mat, const std::vector<std::vector<T> > &src_mat)
{
    assert(dst_mat.size()==src_mat.size());

    for(std::size_t r=0; r<dst_mat.size(); r++)
    {
        assert(dst_mat[r].size()==src_mat[r].size());

        for(std::size_t c=0; c<dst_mat.size(); c++)
        {
            dst_mat[r][c]+=src_mat[r][c];
        }
    }
}

template<typename T>
void mat_sub(std::vector<std::vector<T> > &dst_mat, const std::vector<std::vector<T> > &src_mat)
{
    assert(dst_mat.size()==src_mat.size());

    for(std::size_t r=0; r<dst_mat.size(); r++)
    {
        assert(dst_mat[r].size()==src_mat[r].size());

        for(std::size_t c=0; c<dst_mat.size(); c++)
        {
            dst_mat[r][c]-=src_mat[r][c];
        }
    }
}

template<typename T>
void mat_dot_mult(std::vector<std::vector<T> > &dst_mat, const T &src)
{
    for(std::size_t r=0; r<dst_mat.size(); r++)
    {
        for(std::size_t c=0; c<dst_mat[r].size(); c++)
        {
            dst_mat[r][c]*=src;
        }
    }
}

template<typename T>
void mat_dot_div(std::vector<std::vector<T> > &dst_mat, const T &src)
{
    for(std::size_t r=0; r<dst_mat.size(); r++)
    {
        for(std::size_t c=0; c<dst_mat[r].size(); c++)
        {
            dst_mat[r][c]/=src;
        }
    }
}

template<typename T>
T mat_fabs_sum(const std::vector<std::vector<T> > &src_mat)
{
    T s=0;

    for(std::size_t i = 0; i<src_mat.size(); i++)
    {
        s+=vec_fabs_sum(src_mat[i]);
    }
    return s;
}

template<typename T>
void output_mat(const std::vector<std::vector<T> > &mat)
{
    for(std::size_t i=0; i<mat.size(); i++)
    {
        output_vec(mat[i]);
    }
}

template<typename T>
void vec_mult_mat(const std::vector<T> &vec, const std::vector<std::vector<T> > &mat, std::vector<T> &result)
{
    result.clear();

    assert(vec.size()==mat.size());

    result.resize(mat[0].size(), 0);

    for(std::size_t r=0; r<vec.size(); r++)
    {
        for(std::size_t c=0; c<mat[r].size(); c++)
        {
            result[c]+=vec[r]*mat[r][c];
        }
    }
}

template<typename T>
T mat_determinant(const std::vector<std::vector<T> > &mat);

template<typename T>
T mat_cofactor(const std::vector<std::vector<T> > &mat, const int &row, const int &col)
{
    // calculate the cofactor of mat_{row, col}
    std::vector<std::vector<T> > sub_mat;

    sub_mat.reserve(mat.size()-1);

    for(std::size_t r=0; r<row; r++)
    {
        std::vector<T> sub_row(mat[r].size()-1, 0);
    
        std::copy(mat[r].begin(), mat[r].begin()+col, sub_row.begin());
        
        std::copy(mat[r].begin()+col+1, mat[r].end(), sub_row.begin()+col);

        sub_mat.push_back(sub_row);
    }
    for(std::size_t r=row+1; r<mat.size(); r++)
    {
        std::vector<T> sub_row(mat[r].size()-1, 0);
        
        std::copy(mat[r].begin(), mat[r].begin()+col, sub_row.begin());
        
        std::copy(mat[r].begin()+col+1, mat[r].end(), sub_row.begin()+col);

        sub_mat.push_back(sub_row);
    }
//    std::cout << "Sub matrix: " << mat[row][col] << std::endl;
//    output_mat(sub_mat);

    if((row+col) % 2==0)
        return mat_determinant(sub_mat);
    else
        return -mat_determinant(sub_mat);
}

template<typename T>
T mat_inverse(const std::vector<std::vector<T> > &mat, std::vector<std::vector<T> > &inv_mat)
{
    // calculate the inverse matrix of mat, and return the determinant of mat
    inv_mat.assign(mat.begin(), mat.end());

    if(mat.size()==1)
    {
        inv_mat[0][0]=1.0/mat[0][0];

        return mat[0][0];
    }

    T det=0;

    for(std::size_t row=0; row<mat.size(); row++)
    {
        for(std::size_t col=0; col<mat[row].size(); col++)
        {
            inv_mat[col][row]=mat_cofactor(mat, row, col);

            if(row==0)
            {
                det+=mat[row][col]*inv_mat[col][row];
            }
        }
    }
    mat_dot_div(inv_mat, det);

    return det;
}

template<typename T>
T mat_determinant(const std::vector<std::vector<T> > &mat)
{
    T d;
    if(mat.size()==1)
    {
        return mat[0][0];
    }
    else if(mat.size()==2)
    {
        return mat[0][0]*mat[1][1]-mat[0][1]*mat[1][0];
    }
    else
    {
        d=0;

        for(std::size_t c=0; c<mat[0].size(); c++)
        {
            d+=mat[0][c]*mat_cofactor(mat, 0, c);
        }
    }
    return d;
}

template<typename T>
void output_dim_mat(const std::vector<std::vector<T> > &mat)
{
    std::cout << "Dimension of matrix: " << mat.size() << " * " << mat[0].size() << std::endl;
}

#endif
