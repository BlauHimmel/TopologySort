#ifndef GRAPH_H
#define GRAPH_H

#include "stdafx.h"
#include <stack>
#include <queue>
#include <list>
#include <iostream>
#define INF 2147483647/2
using namespace std;

template <typename VertexType, typename EdgeType>
class Edge
{
public:
	int posOfDestination;	//the position of the another vertex on the edge
	EdgeType Weight;	//weight of the edge
	Edge<VertexType, EdgeType>* link;
	Edge();
	Edge(int dest, EdgeType weight);
	bool operator != (Edge<VertexType, EdgeType>& j) const;
};

template <typename VertexType, typename EdgeType>
class Vertex
{
public:
	VertexType data;
	int inDegree;
	Edge<VertexType, EdgeType>* adj;	//the head pointer of edge chain node
};

template <typename VertexType, typename EdgeType>
class GraphLink				//if the type of return value is integer,which indicates that a position is returned.If VertexType,which indicates that a data stored in vertex is returned.
{
public:
	GraphLink(GraphLink& j);
	GraphLink(int size = 30);
	~GraphLink();
	void Clear();
	VertexType getValue(int position);						//return the value of vertex at the position of "position"
	EdgeType getWeight(int vertex1, int vertex2);			//return the weight of edge between vertex1 and vertex2
	bool insertVertex(VertexType vertex);
	bool removeVertex(int vertex);
	bool insertEdge(int vertex1, int vertex2, EdgeType weight);
	bool removeEdge(int vertex1, int vertex2);
	int getFirstNeighbor(int vertex);						//return the index of fisrt adjacent vertex of "vertex"
	int getNextNeighbor(int vertex, int nextVertex);		//return the index of next adjacent vertext of "nextVertex" which is a adjacent vertex of "vertex"
	bool DFS(GraphLink<VertexType, EdgeType>& graph, const VertexType& vertex, void(*visit)(VertexType data));
	bool BFS(GraphLink<VertexType, EdgeType>& graph, const VertexType& vertex, void(*visit)(VertexType data));
	int getNumberOfVertices();
	int getNumberOfEdges();
	int getVertexPosition(VertexType vertex);		//return the position(index)of the vertex in the graph
	int getInDegree(int position);
	bool TopSort(int* a);
	void ZeroIndegreeEnqueue(queue<VertexType>& Q);
private:
	Vertex<VertexType, EdgeType>* VertexTable;			//array to store vertex
	int maxVertices;		//the maximum number f vertices
	int numEdges;			//the present number of edges
	int numVertices;		//the present number of vertices
};

template <typename VertexType, typename EdgeType>
bool TopologicalSort(GraphLink<VertexType, EdgeType>& Graph, list<int>& node, ofstream& fout);

template <typename VertexType, typename EdgeType>
GraphLink<VertexType, EdgeType>::GraphLink(GraphLink& j)
{
	numEdges = j.numEdges;
	numVertices = j.numVertices;
	maxVertices = j.maxVertices;
	VertexTable = new Vertex<VertexType, EdgeType>[maxVertices];
	for (int i = 0; i < maxVertices; i++)
	{
		VertexTable[i].data = j.getValue(i);
		VertexTable[i].adj = NULL;
		VertexTable[i].inDegree = j.getInDegree(i);
		if (j.VertexTable[i].adj != NULL)
		{
			Edge<VertexType, EdgeType>* pFirst = NULL;

			Edge<VertexType, EdgeType>* pCurrent = j.VertexTable[i].adj;
			Edge<VertexType, EdgeType>* pTemp1 = new Edge<VertexType, EdgeType>;
			pFirst = pTemp1;
			pTemp1->posOfDestination = pCurrent->posOfDestination;
			pTemp1->Weight = pCurrent->Weight;

			Edge<VertexType, EdgeType>* pCurrentThis = pTemp1;
			pCurrent = pCurrent->link;
			
			while (pCurrent != NULL)
			{
				Edge<VertexType, EdgeType>* pTemp = new Edge<VertexType, EdgeType>;
				pTemp->posOfDestination = pCurrent->posOfDestination;
				pTemp->Weight = pCurrent->Weight;

				pCurrentThis->link = pTemp;
				pCurrentThis = pCurrentThis->link;
				pCurrent = pCurrent->link;
			}
			VertexTable[i].adj = pFirst;
		}
	}
}

template <typename VertexType, typename EdgeType>
int GraphLink<VertexType, EdgeType>::getNumberOfEdges()
{
	return numEdges;
}

template <typename VertexType, typename EdgeType>
int GraphLink<VertexType, EdgeType>::getNumberOfVertices()
{
	return numVertices;
}

template <typename VertexType, typename EdgeType>
bool GraphLink<VertexType, EdgeType>::BFS(GraphLink<VertexType, EdgeType>& graph, const VertexType& vertex, void(*visit)(VertexType data))
{
	bool *visited = new bool[numVertices];
	for (int i = 0; i < numVertices; i++)
		visited[i] = false;
	int location = graph.getVertexPosition(vertex);
	if (location == -1)
		return false;

	VertexType temp;
	temp = getValue(location);
	visit(temp);
	visited[location] = true;

	queue<int> myQueue;
	myQueue.push(location);

	int next;
	while (!myQueue.empty())
	{
		location = myQueue.front();
		myQueue.pop();
		next = graph.getFirstNeighbor(location);
		while (next != -1)
		{
			if (visited[next] == false)
			{
				VertexType temp;
				temp = getValue(next);
				visit(temp);
				visited[next] = true;
				myQueue.push(next);
			}
			next = graph.getNextNeighbor(location, next);
		}
	}
	delete[] visited;
	return true;
};

template <typename VertexType, typename EdgeType>
bool GraphLink<VertexType, EdgeType>::DFS(GraphLink<VertexType, EdgeType>& graph, const VertexType& vertex, void(*visit)(VertexType data))
{
	bool *visited = new bool[numVertices];
	for (int i = 0; i < numVertices; i++)
		visited[i] = false;
	int location = graph.getVertexPosition(vertex);
	if (location == -1)
		return false;

	stack<int> myStack;
	myStack.push(location);
	int nextV;

	while (!myStack.empty())
	{
		VertexType temp;
		temp = getValue(myStack.top());

		if (visited[location] == false)
		{
			visit(temp);
			visited[location] = true;
			nextV = graph.getFirstNeighbor(location);
		}
		else
			nextV = graph.getNextNeighbor(location, nextV);

		if (nextV == -1)
		{
			nextV = myStack.top();
			myStack.pop();
			if (!myStack.empty())
				location = myStack.top();
			else
				return true;
		}
		else if (visited[nextV] == false)
		{
			myStack.push(nextV);
			location = nextV;
		}
		else
		{
			nextV = graph.getNextNeighbor(location, nextV);
			if (nextV == -1)
			{
				nextV = myStack.top();
				myStack.pop();
				if (!myStack.empty())
					location = myStack.top();
				else
					return true;
				continue;
			}
			while (visited[nextV] == true)
			{
				nextV = graph.getNextNeighbor(location, nextV);
				if (nextV == -1)
				{
					nextV = myStack.top();
					myStack.pop();
					if (!myStack.empty())
						location = myStack.top();
					else
						return true;
					break;
				}
			}
			if (visited[nextV] == false)
			{
				myStack.push(nextV);
				location = nextV;
			}
		}
	}
	delete[] visited;
	return true;
};

template <typename VertexType, typename EdgeType>
int GraphLink<VertexType, EdgeType>::getVertexPosition(VertexType vertex)
{
	for (int i = 0; i < numVertices; i++)
	{
		if (VertexTable[i].data == vertex)
			return i;
	}
	return -1;
};


template <typename VertexType, typename EdgeType>
bool GraphLink<VertexType, EdgeType>::removeEdge(int vertex1, int vertex2)
{
	if (vertex1 != -1 && vertex2 != -1)		//legal parameter
	{
		Edge<VertexType, EdgeType> *current = VertexTable[vertex1].adj, *previous = NULL;
		while (current != NULL&&current->posOfDestination != vertex2)	//find the edge to be deleted in vertex1
		{
			previous = current;
			current = current->link;
		}
		if (current != NULL)		//if found
		{
			if (current == VertexTable[vertex1].adj)			//if it is the first adjacent vertex
				VertexTable[vertex1].adj = current->link;
			else
				previous->link = current->link;
			delete current;
			VertexTable[vertex2].inDegree--;
		}
		else	//target not exist
		{
			return false;
		}
		return true;
	}
	return false;
};

template <typename VertexType, typename EdgeType>
bool GraphLink<VertexType, EdgeType>::insertEdge(int vertex1, int vertex2, EdgeType weight)
{
	if (vertex1 >= 0 && vertex2 >= 0 && vertex1 < numVertices&&vertex2 < numVertices)	//legal parameter
	{
		Edge<VertexType, EdgeType> *q = NULL, *p = VertexTable[vertex1].adj;
		while (p != NULL&&p->posOfDestination != vertex2)	//search adjacent vertex "vertex2"
			p = p->link;
		if (p != NULL)
			return false;		//an edge existed between vertex1 and vertex2
		p = new Edge<VertexType, EdgeType>;

		p->posOfDestination = vertex2;
		p->link = VertexTable[vertex1].adj;
		p->Weight = weight;
		VertexTable[vertex1].adj = p;

		VertexTable[vertex2].inDegree++;

		numEdges++;
		return true;
	}
	return false;
};

template <typename VertexType, typename EdgeType>
bool GraphLink<VertexType, EdgeType>::removeVertex(int vertex)
{
	if (numVertices == 0 || vertex == -1 || vertex >= numVertices)	//not enought vertices,or the vertex is out of range
		return false;

	Edge<VertexType, EdgeType> *del = NULL, *pre = NULL, *cur = NULL;

	while (VertexTable[vertex].adj != NULL)
	{
		del = VertexTable[vertex].adj;
		cur = del;
		VertexTable[vertex].adj = cur->link;
		VertexTable[del->posOfDestination].inDegree--;
		delete del;
		numEdges--;
	}

	/*ERROR*/
	for (int i = 0; i < numVertices; i++)
	{
		pre = NULL;
		del = VertexTable[i].adj;
		cur = del;
		while (del != NULL&&del->posOfDestination != vertex)
		{
			pre = del;
			del = del->link;
			cur = del;
		}
		if (cur == NULL)
			continue;
		if (pre == NULL)
		{
			VertexTable[i].adj = cur->link;
			delete del;
		}
		else
		{
			pre->link = cur->link;
			delete del;
		}
	}

	Edge<VertexType, EdgeType> *present = NULL, *target = NULL;
	numVertices--;
	VertexTable[vertex].data = VertexTable[numVertices].data;
	VertexTable[vertex].adj = VertexTable[numVertices].adj;
	VertexTable[vertex].inDegree = VertexTable[numVertices].inDegree;
	for (int i = 0; i < numVertices; i++)
	{
		present = VertexTable[i].adj;
		while (present != NULL&&present->posOfDestination != numVertices)
		{
			present = present->link;
		}
		if (present == NULL)
			continue;
		if (present->posOfDestination == numVertices)
		{
			present->posOfDestination = vertex;
		}

	}	
	return true;
};

template <typename VertexType, typename EdgeType>
bool GraphLink<VertexType, EdgeType>::insertVertex(VertexType vertex)
{
	if (numVertices == maxVertices)
		return false;
	VertexTable[numVertices].data = vertex;
	numVertices++;
	return true;
};

template <typename VertexType, typename EdgeType>
EdgeType GraphLink<VertexType, EdgeType>::getWeight(int vertex1, int vertex2)
{
	if (vertex1 != -1 && vertex2 != -1)
	{
		Edge<VertexType, EdgeType>* temp = VertexTable[vertex1].adj;
		while (temp != NULL&&temp->posOfDestination != vertex2)
			temp = temp->link;
		if (temp != NULL)
			return temp->Weight;
	}
	return INF;
};

template <typename VertexType, typename EdgeType>
int GraphLink<VertexType, EdgeType>::getInDegree(int position)
{
	if (position >= 0 && position < numVertices)
		return VertexTable[position].inDegree;
}

template <typename VertexType, typename EdgeType>
VertexType GraphLink<VertexType, EdgeType>::getValue(int position)
{
	if (position >= 0 && position < numVertices)
		return VertexTable[position].data;
//	else
//		return NULL;
};

template <typename VertexType, typename EdgeType>
int GraphLink<VertexType, EdgeType>::getNextNeighbor(int vertex, int nextVertex)
{
	if (vertex != -1)
	{
		Edge<VertexType, EdgeType>* temp = VertexTable[vertex].adj;
		while (temp != NULL&&temp->posOfDestination != nextVertex)
			temp = temp->link;
		if (temp != NULL&&temp->link != NULL)
			return temp->link->posOfDestination;
	}
	return -1;		//not existed!
};

template <typename VertexType, typename EdgeType>
int GraphLink<VertexType, EdgeType>::getFirstNeighbor(int vertex)
{
	if (vertex != -1)		//if existed!
	{
		Edge<VertexType, EdgeType>* temp = VertexTable[vertex].adj;
		if (temp != NULL)
			return temp->posOfDestination;
	}
	return -1;	//not existed
};

template <typename VertexType, typename EdgeType>
GraphLink<VertexType, EdgeType>::~GraphLink()
{
	for (int i = 0; i < numVertices; i++)
	{
		Edge<VertexType, EdgeType> *temp = VertexTable[i].adj;
		while (temp != NULL)
		{
			VertexTable[i].adj = temp->link;
			delete temp;
			temp = VertexTable[i].adj;
		}
	}
	delete[] VertexTable;
};

template <typename VertexType, typename EdgeType>
void GraphLink<VertexType, EdgeType>::Clear()
{
	for (int i = 0; i < numVertices; i++)
	{
		Edge<VertexType, EdgeType> *temp = VertexTable[i].adj;
		while (temp != NULL)
		{
			VertexTable[i].adj = temp->link;
			delete temp;
			temp = VertexTable[i].adj;
		}
		VertexTable[i].inDegree = 0;
	}
	numVertices = 0;
	numEdges = 0;
};

template <typename VertexType, typename EdgeType>
GraphLink<VertexType, EdgeType>::GraphLink(int size)
{
	maxVertices = size;
	numVertices = 0;
	numEdges = 0;
	VertexTable = new Vertex<VertexType, EdgeType>[maxVertices];
	for (int i = 0; i < maxVertices; i++)
	{
		VertexTable[i].adj = NULL;
		VertexTable[i].inDegree = 0;
	}
};

template <typename VertexType, typename EdgeType>
bool Edge<VertexType, EdgeType>::operator != (Edge<VertexType, EdgeType>& j) const
{
	return (posOfDestination != j.posOfDestination) ? true : false;
};

template <typename VertexType, typename EdgeType>
Edge<VertexType, EdgeType>::Edge()
{
	link = NULL;
};

template <typename VertexType, typename EdgeType>
Edge<VertexType, EdgeType>::Edge(int dest, EdgeType weitht)
{
	posOfDestination = dest;
	Weight = weight;
	link = NULL;
};

template <typename VertexType, typename EdgeType>
bool GraphLink<VertexType, EdgeType>::TopSort(int* a)
{
	for (int i = 0; numVertices > 0; i++)
	{
		int zero;
		for (zero = 0; zero < numVertices; zero++)
		{
			if (VertexTable[zero].inDegree == 0)
				break;
		}
		if (numVertices != 0 && zero >= numVertices)
		{
			return false;
		}
		a[i] = VertexTable[zero].data;
		removeVertex(zero);
	}
	return true;
}

template <typename VertexType, typename EdgeType>
void GraphLink<VertexType, EdgeType>::ZeroIndegreeEnqueue(queue<VertexType>& Q)
{
	for (int i = 0; i < numVertices; i++)
	{
		if (VertexTable[i].inDegree == 0)
		{
			Q.push(VertexTable[i].data);
		}
	}
}

template <typename VertexType, typename EdgeType>
bool TopologicalSort(GraphLink<VertexType, EdgeType>& Graph, list<int>& node, ofstream& fout)
{
	static bool bFinish = false;
	queue<int> Q;
	Graph.ZeroIndegreeEnqueue(Q);
	while (!Q.empty()) 
	{
		int nVertex = Q.front();
		Q.pop();
		list<int> copy_node;

		copy_node.assign(node.begin(), node.end());
		copy_node.push_back(nVertex);

		GraphLink<VertexType, EdgeType> CGraphCopy(Graph);
		CGraphCopy.removeVertex(CGraphCopy.getVertexPosition(nVertex));

		if (CGraphCopy.getNumberOfVertices() == 0)
		{
			while (!copy_node.empty())
			{
				int nTemp = copy_node.front();
				copy_node.pop_front();
				fout << "P" << nTemp;
				if (copy_node.size() == 0)
					break;
				fout << ",";
			}
			fout << endl;
			bFinish = true;
		}
		else {
			TopologicalSort(CGraphCopy, copy_node, fout);
		}
	}
	return bFinish;
}

#endif
