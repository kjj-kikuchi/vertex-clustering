
// QEMを用いた頂点クラスタリングによるメッシュ簡略化

#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <numbers>
#include <chrono>
#include <Eigen/Dense>

struct Mesh
{
    std::vector<Eigen::Vector3d> V;
    std::vector<Eigen::Vector3i> F;
};

struct EigenVector3iHasher
{
    std::size_t operator()(Eigen::Vector3i const& key) const
    {
        std::hash<int> hasher;
        std::size_t seed = 0;
        seed ^= hasher(key(0)) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= hasher(key(1)) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= hasher(key(2)) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

class VertexClustering
{
    std::vector<Eigen::Vector3d> const *V;
    std::vector<Eigen::Vector3i> const *F;

    Eigen::Vector3d min;            // メッシュの最小のx,y,z座標
    Eigen::Vector3d max;            // メッシュの最大のx,y,z座標
    double longest_edge_length;     // bbの最長辺の長さ
    double cell_length;             // セルの一辺の長さ

    struct QuadricErrorMetric
    {
        Eigen::Matrix3d A;
        Eigen::Vector3d b;

        Eigen::Vector3d compute(double weight, Eigen::Vector3d c)
        {
            Eigen::Matrix3d A_dash = A + weight * Eigen::Matrix3d::Identity();
            Eigen::Vector3d b_dash = b + weight * c;
            return A_dash.inverse() * b_dash;
        }
    };

public:
    std::vector<Eigen::Vector3d> newV;
    std::vector<Eigen::Vector3i> newF;

    Eigen::Vector3i cell_num;       // セルの個数 (nx, ny, nz)


    VertexClustering(Mesh *mesh, double cell_length_percentage)
    {
        V = &(mesh->V);
        F = &(mesh->F);
        cell_length = cell_length_percentage;
    }

    void calc_bounding_box()
    {
        min = (*V)[0];
        max = (*V)[0];

        for (int i = 0; i < (*V).size(); i++)
        {
            for (int j = 0; j < 3; j++)
            {
                if ((*V)[i](j) < min(j)) min(j) = (*V)[i](j);
                if ((*V)[i](j) > max(j)) max(j) = (*V)[i](j);
            }
        }
        longest_edge_length = std::max( {max(0) - min(0), max(1) - min(1), max(2) - min(2)} );
    }

    void get_cell_number()
    {
        cell_length *= longest_edge_length;
        for (int i = 0; i < 3; i++)
        {
            cell_num(i) = std::ceil((max(i) - min(i)) / cell_length) + 1;
        }
    }

    void compute()
    {
        get_cell_number();  // セルの個数 nx, ny, nz を取得

        // 面リストを走査し，頂点のセル番号を調べる
        // 3頂点のセル番号が全て異なるとき，マップに各セルを追加・面リストを追加する
        std::unordered_map<Eigen::Vector3i, int, EigenVector3iHasher> cell_map; // Key : セル番号(ix iy iz), Value : 新しい頂点番号
        std::vector<int> cell_v_num;                                            // セルに含まれる頂点の個数
        int newV_idx = 0;                                                       // 新しい頂点番号

        for (auto& f : (*F))
        {
            std::vector<Eigen::Vector3i> cell_idx(3);
            cell_idx[0] = ( ((*V)[f(0)] - min) / cell_length ).cast<int>();
            cell_idx[1] = ( ((*V)[f(1)] - min) / cell_length ).cast<int>();
            cell_idx[2] = ( ((*V)[f(2)] - min) / cell_length ).cast<int>();

            if (cell_idx[0] != cell_idx[1] &&
                cell_idx[1] != cell_idx[2] &&
                cell_idx[2] != cell_idx[0])
            {
                Eigen::Vector3i face;
                for (int i = 0; i < 3; i++)
                {
                    if (! cell_map.contains(cell_idx[i]))
                    {
                        cell_map.emplace(cell_idx[i], newV_idx);
                        newV_idx ++;
                    }

                    face(i) = cell_map.at(cell_idx[i]);
                }
                newF.push_back(face);
            }
        }

        // セル内の頂点集合の重心を計算
        // std::vector<Eigen::Vector3d> cell_centroid(cell_map.size(), Eigen::Vector3d::Zero());   // セル内の重心
        cell_v_num.resize(cell_map.size(), 0);
        for (int i = 0; i < (*V).size(); i++)
        {
            Eigen::Vector3i cell_idx = ( ((*V)[i] - min) / cell_length ).cast<int>();
            if (cell_map.contains(cell_idx))
            {
                //cell_centroid[ cell_map.at(cell_idx) ] += (*V)[i];
                cell_v_num[ cell_map.at(cell_idx) ] ++;
            }
        }
//        for (auto& [cell, idx] : cell_map)
//        {
//            cell_centroid[idx] /= cell_v_num[idx];
//        }

        //
        std::vector<QuadricErrorMetric> QEM(cell_map.size());
        for (auto& f : (*F))
        {
            Eigen::Vector3d normal = ( ( (*V)[f(1)] - (*V)[f(0)] ).cross( (*V)[f(2)] - (*V)[f(0)] ) ).normalized();
            for (int i = 0; i < 3; i++)
            {
                Eigen::Vector3i cell_idx = ( ((*V)[f(i)] - min) / cell_length ).cast<int>();
                if (cell_map.contains(cell_idx))
                {
                    int idx = cell_map.at(cell_idx);
                    Eigen::Matrix3d A = normal * normal.transpose();
                    QEM[idx].A += A;
                    QEM[idx].b += A * (*V)[f(i)];
                }
            }
        }

        // 各セルのQEMを計算し，新しい頂点座標を計算
        newV.resize(cell_map.size(), Eigen::Vector3d::Zero());
        for (auto& [cell, idx] : cell_map)
        {
            Eigen::Vector3d c = cell.cast<double>() * cell_length + min + Eigen::Vector3d(0.5, 0.5, 0.5) * cell_length;
            double weight = cell_v_num[idx] * 0.001;
//            newV[idx] = QEM[idx].compute(weight, cell_centroid[idx]);
            newV[idx] = QEM[idx].compute(weight, c);
        }
    }
};



void read_obj(std::string const& filename, Mesh& mesh)
{
    std::ifstream ifs(filename);
    if (ifs.fail())
    {
        std::cerr << "Failed to open file." << "\n";
        std::exit(1);
    }

    std::string line;
    while (std::getline(ifs, line))
    {
        std::istringstream string_in(line);
        std::string type;
        string_in >> type;

        if (type == "v")
        {
            Eigen::Vector3d v;
            string_in >> v(0) >> v(1) >> v(2);
            mesh.V.push_back(v);
        }
        else if (type == "f")
        {
            Eigen::Vector3i f;
            string_in >> f(0) >> f(1) >> f(2);
            f -= Eigen::Vector3i{1, 1, 1};
            mesh.F.push_back(f);
        }
    }
}

void export_obj(std::string name, std::vector<Eigen::Vector3d> const& vert, std::vector<Eigen::Vector3i> const& face)
{
    std::ofstream of;
    name.erase(name.length()-4);
    std::string filename = name + std::to_string(vert.size()) + "_simplified.obj";
    of.open(filename, std::ios::out);
    for(auto& v : vert)
    {
        of << "v " << v(0) << " " << v(1) << " " << v(2) << std::endl;
    }
    for(auto& f : face)
    {
        of << "f " << f(0)+1 << " " << f(1)+1 << " " << f(2)+1 << std::endl;
    }
    of.close();
}






int main(int argc, char *argv[])
{
    // Enter cell length percentage .........................
    double d;
    std::cout << "Enter percentage of cell length : ";
    std::cin >> d;

    // Read mesh file .......................................
    auto start_readmesh = std::chrono::system_clock::now();

    std::string filename;
    if (argc != 2)
    {
        std::cout << "wrong command line argument" << std::endl;
        std::exit(1);
    }
    filename = std::string(argv[1]);

    Mesh mesh;
    read_obj(filename, mesh);

    auto end_readmesh = std::chrono::system_clock::now();


    // Vertex Clustering ....................................
    auto start_vc = std::chrono::system_clock::now();

    VertexClustering vc(&mesh, d);
    vc.calc_bounding_box();
    vc.compute();

    auto end_vc = std::chrono::system_clock::now();

    // Export and output .....................................
    export_obj(filename, vc.newV, vc.newF);

    using namespace std::chrono_literals;
    std::cout << "cell number : " << vc.cell_num.transpose() << ", ";
    std::cout << vc.cell_num(0) * vc.cell_num(1) * vc.cell_num(2) << "\n\n";
    std::cerr << "Read mesh : " << (end_readmesh - start_readmesh) / 1.0s << " sec\n";
    std::cerr << "Vertex clustering : " << (end_vc - start_vc) / 1.0s << " sec\n";
    //std::cerr << "Export mesh : " << (end_export - start_export) / 1.0s << " sec\n";

    return 0;
}

