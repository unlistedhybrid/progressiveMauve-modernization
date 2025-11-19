#include "AlignmentTree.h"
#include <stack>
#include <sstream>
using namespace std;

AlignmentTree::AlignmentTree() : vector<TreeNode>() {
	root = 0;
}

AlignmentTree::~AlignmentTree() {}

void AlignmentTree::clear() {
	vector<TreeNode>::clear();
	root = 0;
}

void AlignmentTree::readTree(istream& tree_file) {
	string line;
	clear();
	if(!getline(tree_file, line))
		return;

	stringstream line_str(line);
	string tree_line;
	getline(line_str, tree_line, ';');
	stack<node_id_t> node_stack;
	stringstream blen_str;
	TreeNode new_node;
	new_node.distance = 0;
	
	// Fixed: uint -> unsigned int
	for(unsigned int charI = 0; charI < tree_line.size(); charI++) {
		switch(tree_line[charI]) {
			case '(':
				if(node_stack.size() > 0) {
					new_node.parents.clear();
					new_node.parents.push_back(node_stack.top());
					(*this)[node_stack.top()].children.push_back((node_id_t)(*this).size());
				}
				node_stack.push((node_id_t)(*this).size());
				push_back(new_node);
				break;
			case ')':
				node_stack.pop();
				break;
			case ':':
				break;
			default:
				break;
		}
	}
}

void AlignmentTree::writeTree(ostream& os) const {
	stack<node_id_t> node_stack;
	stack<unsigned int> child_stack; // Fixed: uint -> unsigned int
	node_stack.push(root);
	child_stack.push(0);
	os << "(";
	while(node_stack.size() > 0) {
		if((*this)[node_stack.top()].children.size() != 0) {
			if(child_stack.top() == (*this)[node_stack.top()].children.size()) {
				os << ")";
				node_stack.pop();
				child_stack.pop();
				continue;
			}
			node_id_t child = (*this)[node_stack.top()].children[child_stack.top()];
			node_stack.push(child);
			child_stack.top()++;
			if(child_stack.top() > 1)
				os << ",";
			if((*this)[child].children.size() > 0)
				child_stack.push(0);
			continue;
		}
		os << (*this)[node_stack.top()].name << ":" << (*this)[node_stack.top()].distance;
		node_stack.pop();
	}
	os << ";" << endl;
}

double AlignmentTree::getHeight() const {
	return getHeight(root);
}

double AlignmentTree::getHeight(node_id_t nodeI) const {
	if((*this)[nodeI].children.size() == 0)
		return (*this)[nodeI].distance;
	return (*this)[nodeI].distance + getHeight((*this)[nodeI].children[0]);
}
