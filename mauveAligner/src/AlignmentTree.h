#ifndef ALIGNMENTTREE_H
#define ALIGNMENTTREE_H

#include <vector>
#include <string>
#include <iostream>

typedef unsigned int node_id_t;

struct TreeNode {
    std::string name;
    double distance;
    std::vector<node_id_t> children;
    std::vector<node_id_t> parents;

    TreeNode() : name(""), distance(0.0) {}
};

class AlignmentTree : public std::vector<TreeNode> {
public:
    AlignmentTree();
    ~AlignmentTree();

    void clear();
    void readTree(std::istream& tree_file);
    void writeTree(std::ostream& os) const;

    double getHeight() const;
    double getHeight(node_id_t nodeI) const;

    node_id_t root;
};

#endif
