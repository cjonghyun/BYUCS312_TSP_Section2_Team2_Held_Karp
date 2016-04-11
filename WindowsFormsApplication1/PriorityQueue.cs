using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TSP
{
    /**
    * Priority Queue class
    * Binary Heap is used to create a priority queue.
    Time complexity for add: average O (log n)
    Time complexity for remove: average O(log n)
    **/    class PriorityQueue
    {
        private List<ReducedMatrix> tree;
        private int size;
        public PriorityQueue()
        {
            tree = new List<ReducedMatrix>();
            tree.Add(null);
            size = 0;
        }
        public int getSize()
        {
            return size;
        }
        // swap two nodes
        public void swap(int a, int b)
        {
            ReducedMatrix temp = tree[a];
            tree[a] = tree[b];
            tree[b] = temp;
        }
        /**
        add a new node.
            Insert the new node as a leaf node (last index)
            and start comparing with its root nodes and swap positions if the new node has the small bound.
        **/
        public void add(ReducedMatrix newMatrix)
        {            
            tree.Add(newMatrix);
            size++;
            int k = size;
            while (k > 1)
            {
                if (tree[k / 2].bound > tree[k].bound)
                {
                    swap(k, k / 2);
                }
                else
                    break;
                k = k / 2;
            }
        }

        /**
        remove function for priority queue
            save the root node(the node with the minimum bound) and put the node with the last index as a root node.
            start comparing the bound values with its children and find the right place for this node.
        **/
        public ReducedMatrix remove()
        {
            if (size ==0)
                return null;
            ReducedMatrix root = tree[1];
            tree[1] = tree[size];
            tree.RemoveAt(size);
            size--;
            int k = 1;
            Boolean flag = true;
            if(size == 0)
            {
                return root;
            }
            while (flag)
            {
                if (k * 2 > size)
                    break;
                flag = false;
                if(k*2 == size)
                {
                    if (tree[k * 2].bound > tree[k].bound)
                    {
                        swap(k, k * 2);
                        k *= 2;
                        flag = true;
                    }
                    break;
                }
                if (tree[k * 2].bound < tree[k * 2 + 1].bound)
                {
                    if (tree[k * 2].bound < tree[k].bound)
                    {
                        swap(k, k * 2);
                        k *= 2;
                        flag = true;
                    }
                    else if (tree[k * 2 + 1].bound < tree[k].bound)
                    {
                        swap(k * 2 + 1, k);
                        k = k * 2 + 1;
                        flag = true;
                    }
                }
                else
                {
                    if (tree[k * 2 + 1].bound < tree[k].bound)
                    {
                        swap(k, k * 2 + 1);
                        k = k * 2 + 1;
                        flag = true;
                    }
                    else if(tree[k*2].bound < tree[k].bound)
                    {
                        swap(k, k * 2);
                        k = k * 2;
                        flag = true;
                    }
                }
            }
            return root;
        }        
    }
}
