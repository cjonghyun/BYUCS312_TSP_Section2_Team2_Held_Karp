using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.Serialization.Formatters.Binary;
using System.Text;
using System.Threading.Tasks;



namespace TSP
{
    /**
    class used to store the information of the state including matrix
        stores the reduced matrix and bound, number of the current city and the route of the tour so far.
        It is used to calculate the bound by getting reduced matrix
    **/
    class ReducedMatrix
    {

        /**
        a method for doing deep copy the objects from the parent state.
            There were some cases that information of states are changed when the child state's information is changed. (for example, reducedMatrix )
            In order to prevent that problem, deep copy is used
        **/
        public static T DeepCopy<T>(T other)
        {
            using (MemoryStream ms = new MemoryStream())
            {
                BinaryFormatter formatter = new BinaryFormatter();
                formatter.Serialize(ms, other);
                ms.Position = 0;
                return (T)formatter.Deserialize(ms);
            }
        }

        public double[,] matrix;
        public double bound;
        public int n;
        public int cityNumber;
        public ArrayList Route;

        /**
        Constructor for Creating the matrix for the start node
        **/
        public ReducedMatrix(double[,] m, int size, City city)
        {
            matrix = m;
            bound = 0;
            n = size;
            Route = new ArrayList();
            cityNumber = 0;
            Route.Add(city);
        }

        
        public ReducedMatrix(ReducedMatrix parent, City nextCity, int citynum)
        {
            matrix = DeepCopy<double[,]>(parent.matrix);
            bound = parent.bound + parent.matrix[parent.cityNumber, citynum];
            Route = new ArrayList(parent.Route);
            cityNumber = citynum;
            n = parent.n;
            Route.Add(nextCity);
        }

        public void cancel(int row, int column)
        {
            for (int i = 0; i < n; i++)
            {
                matrix[row, i] = double.PositiveInfinity;
                matrix[i, column] = double.PositiveInfinity;
            }
        }

        /**
        Method for getting ReducedMatrix

        **/
        public void reduce()
        {
            Boolean flag;
            double min;
            int k;
            for (int i = 0; i < n; i++)
            {
                min = double.PositiveInfinity;
                k = -1;
                for (int j = 0; j < n; j++)
                {
                    if (matrix[i, j] < min)
                    {
                        min = matrix[i, j];
                        k = j;
                        if (min == 0)
                            break;
                    }
                }
                if (min != 0 && min != double.PositiveInfinity)
                {
                    for (int j = 0; j < n; j++)
                    {
                        if (matrix[i, j] != double.PositiveInfinity && matrix[i,j] != 0)
                            matrix[i, j] -= min;
                    }
                    bound += min;
                }
            }

            for (int i = 0; i < n; i++)
            {
                flag = false;
                for (int j = 0; j < n; j++)
                {
                    if (matrix[j, i] == 0)
                        flag = true;
                }
                if (!flag)
                {
                    min = double.PositiveInfinity;
                    k = -1;
                    for (int j = 0; j < n; j++)
                    {
                        if (matrix[j, i] < min)
                        {
                            min = matrix[j, i];
                            k = j;
                            if (min == 0)
                                break;
                        }
                    }
                    if (min!=0 && min != double.PositiveInfinity)
                    {
                        for (int j = 0; j < n; j++)
                        {
                            if (matrix[j, i] != double.PositiveInfinity && matrix[j,i] != 0)
                            {
                                matrix[j, i] -= min;
                            }
                        }
                        bound += min;
                    }
                }
            }
        }
    }
}