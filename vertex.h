#ifndef VERTEX_H
#define VERTEX_H
#include <iostream>
#include <memory>
#include <vector>
#include <unordered_map>
#include <stdio.h>
#include <assert.h>
#include "Eigen/Dense"

#define V(id) \
vertices[id] = std::make_shared<vertex>(id);

struct vertex{
    size_t id; //id in the subdivided graph
    std::vector<size_t> adj; //adjacent table
    vertex()=delete;
    vertex(size_t id):id(id){}
    void print(){
        std::vector<size_t> tmp = adj;
        sort(tmp.begin(), tmp.end(), std::less<size_t>());
        printf("%zu: ", id);
        for(size_t v: tmp){
            printf("%zu ", v);
        }
        printf("\n");
    }
};

struct UnaryFct {
   double operator()(float x) const { std::cout << "HERE\n"; return x < newthres?0.0f:1.0f; }
   mutable float newthres;
};

struct helper{
    static std::shared_ptr<Eigen::MatrixXf> construct(Eigen::MatrixXf& original, const size_t k, \
        std::unordered_map<size_t, std::shared_ptr<vertex>>& vertices){
        //original: the Laplacian of the original graph
        //Output the Laplacian of the subdivided graph

        //Note: Integer type `Eigen::MatrixXi` works well for math but poorly here, 
        //as some important functions, e.g., determinant(), do not support the type `Eigen::MatrixXi`.
        

        //Clear the vertices first
        vertices.clear();
        
        //Check:
        if(k <= 1 || original.cols() != original.rows()){
            return nullptr;
        }
        size_t id = original.cols();
        for(size_t j = 0; j < original.cols(); ++j){
            V(j); //new
        }

        //Construct adjacent table
        for(size_t j = 0; j < original.cols(); ++j){//Eigen is column-first 
            for(size_t i = j+1; i < original.rows(); ++i){
                if(original(i, j) == -1){ //create new subdivision
                    if(original(j, i) != -1){
                        return nullptr; //Check symmetry
                    }
                    if(k == 2){
                        V(id);
                        vertices[i] -> adj.push_back(id);
                        vertices[j] -> adj.push_back(id);
                        vertices[id] -> adj.push_back(i);
                        vertices[id] -> adj.push_back(j);
                        ++id;
                    }
                    else{
                        for(size_t idtmp = id; idtmp <= id + k - 2; ++idtmp){
                            V(idtmp);
                            if(idtmp == id){
                                //connect idtmp to j, which is minimum between {i, j}
                                vertices[j] -> adj.push_back(idtmp);
                                vertices[idtmp] -> adj.push_back(j);
                                vertices[idtmp] -> adj.push_back(k == 2 ? i : idtmp+1);
                            }
                            else if(idtmp == id + k - 2){
                                //connect idtmp to i, which is maximum between {i, j}
                                vertices[i] -> adj.push_back(idtmp);
                                vertices[idtmp] -> adj.push_back(i);
                                vertices[idtmp] -> adj.push_back(k == 2 ? j : idtmp - 1);
                            }
                            else{
                                vertices[idtmp] -> adj.push_back(idtmp - 1);
                                vertices[idtmp] -> adj.push_back(idtmp + 1);
                            }
                        }
                        id += (k - 1);
                    }
                }
            }//for
        }

        //Check diagonal element:
        for(size_t j = 0; j < original.cols(); ++j){
            if(!vertices[j] || vertices[j] -> adj.size() != original(j, j)){
                return nullptr;
            }
        }

        //Construct Laplacian matrix (of the subdivided graph)
        assert(id == vertices.size());
        std::shared_ptr<Eigen::MatrixXf> result(new Eigen::MatrixXf(id, id));//Why does `Eigen::MatrixXf::Zero` not work???
        result -> setZero();
        for(size_t i = 0; i < id; ++i){
            (*result)(i, i) = (int)vertices[i] -> adj.size();
            for(size_t v: vertices[i] -> adj){
                (*result)(i, v) = (*result)(v, i) = -1.0f;
            }
        }

        return result;
    }

    static std::shared_ptr<Eigen::MatrixXf> construct(std::shared_ptr<Eigen::MatrixXf> original, const size_t k, \
        std::unordered_map<size_t, std::shared_ptr<vertex>>& vertices){
            return (original ? construct(*original, k, vertices) : nullptr);
        }

    

    void print(std::unordered_map<size_t, std::shared_ptr<vertex>>& vertices){
        std::cout << "Totally " << vertices.size() << " vertices.\n";
        for(auto it = vertices.begin(); it != vertices.end(); ++it){
            if(it -> second){
                it -> second -> print();
            }
        }
    }

    static size_t Kirchhoff(Eigen::MatrixXf Laplacian){
        //Kirchhoff Matrix's tree theorem: https://en.wikipedia.org/wiki/Kirchhoff%27s_theorem
        //subtract arbitrary minor
        if(Laplacian.rows() != Laplacian.cols() || Laplacian.rows() < 1){
            return -1;
        }
        Eigen::VectorXi indices(Laplacian.rows()-1);
        for(int i = 0; i <= Laplacian.rows()-2; ++i){
            indices(i) = i;
        }
        Eigen::MatrixXf newmatrix = Laplacian(indices, indices);
        return (size_t)newmatrix.determinant();
    }

    static size_t Kirchhoff(std::shared_ptr<Eigen::MatrixXf> Laplacian){
        if(!Laplacian){
            return -1;
        }
        return Kirchhoff(*Laplacian);
    }

    static std::shared_ptr<Eigen::MatrixXf> genKn(size_t n){
        if(n <= 1){
            return nullptr;
        }
        std::shared_ptr<Eigen::MatrixXf> lap(new Eigen::MatrixXf(n, n)); 
        for(size_t i = 0; i < n; ++i){
            for(size_t j = 0; j < n; ++j){
                (*lap)(i, j) = ((i == j) ? (n - 1) : -1.0f);
            }
        }
        return lap;
    }

    static std::shared_ptr<Eigen::MatrixXf> genKmn(size_t m, size_t n){
        //Biparite
        if(m < 1 || n < 1){
            return nullptr;
        }
        std::shared_ptr<Eigen::MatrixXf> lap(new Eigen::MatrixXf(m+n, m+n)); 
        for(size_t i = 0; i < m+n; ++i){
            for(size_t j = 0; j < m+n; ++j){
                if(i == j){
                    if(i < m){
                        (*lap)(i, j) = n;
                    }else{
                        (*lap)(i, j) = m;
                    }
                }else{
                    if(i < m && j >= m){
                        (*lap)(i, j) = -1.0f;
                    }
                    else if(i >= m && j < m){
                        (*lap)(i, j) = -1.0f;
                    }
                    else{
                        (*lap)(i, j) = 0.0f;
                    }
                }
            }
        }
        return lap;
    }

    static std::shared_ptr<Eigen::MatrixXf> genrandom(size_t n, float thres=0.5){
        if(n < 2){
            return nullptr;
        }
        std::shared_ptr<Eigen::MatrixXf> lap(new Eigen::MatrixXf(n, n)); 
        Eigen::MatrixXf tmpr(n, n);
        Eigen::MatrixXf tmp(n, n);
        srand((unsigned int) time(0));
        tmpr.setRandom();
        tmp.setConstant(2.0f);
        *lap = 0.25*(tmpr + tmpr.transpose() + tmp);
        for(size_t i = 0; i < n; ++i){
            float sumi = 0.0f;
            for(size_t j = 0; j < n; ++j){
                if(i != j){
                    (*lap)(i, j) = ((*lap)(i, j) < thres) ? 0.0f : -1.0f;
                    sumi += (*lap)(i, j);
                }
            }
            (*lap)(i, i) = -sumi;
        }
        return lap;
    }
};

#endif