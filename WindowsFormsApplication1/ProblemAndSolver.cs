using System;
using System.Collections;
using System.Collections.Generic;
using System.Text;
using System.Drawing;
using System.Diagnostics;
using System.Windows.Forms;
using System.Linq;

namespace TSP
{

    class ProblemAndSolver
    {

        private class TSPSolution
        {
            /// <summary>
            /// we use the representation [cityB,cityA,cityC] 
            /// to mean that cityB is the first city in the solution, cityA is the second, cityC is the third 
            /// and the edge from cityC to cityB is the final edge in the path.  
            /// You are, of course, free to use a different representation if it would be more convenient or efficient 
            /// for your data structure(s) and search algorithm. 
            /// </summary>
            public ArrayList
                Route;
            /// <summary>
            /// constructor
            /// </summary>
            /// <param name="iroute">a (hopefully) valid tour</param>
            public TSPSolution(ArrayList iroute)
            {
                Route = new ArrayList(iroute);
            }

            /// <summary>
            /// Compute the cost of the current route.  
            /// Note: This does not check that the route is complete.
            /// It assumes that the route passes from the last city back to the first city. 
            /// </summary>
            /// <returns></returns>
            public double costOfRoute()
            {
                // go through each edge in the route and add up the cost. 
                int x;
                City here;
                double cost = 0D;

                for (x = 0; x < Route.Count - 1; x++)
                {
                    here = Route[x] as City;
                    cost += here.costToGetTo(Route[x + 1] as City);
                }

                // go from the last city to the first. 
                here = Route[Route.Count - 1] as City;
                cost += here.costToGetTo(Route[0] as City);
                return cost;
            }
        }
     

        private class Node {
               public int Vertex { get; set; }
               public Node[] ChildNodes { get; set; }
               public bool Selected { get; set; }
        }

        public class MNode :IEquatable<MNode> 
        {
            public int Vertex { get; set; }
            public HashSet<int> citySet { get; set; }

            public override int GetHashCode() {
                return Vertex * 100000 + citySet.Sum();
            }

            public override bool Equals(Object other)
            {
                return Equals(other as MNode);
            }
            public  bool Equals(MNode other)
            {
                if (Vertex != other.Vertex)
                    return false;
                if(Vertex == other.Vertex && citySet.Count == other.citySet.Count)
                {
                    int[] tempA = citySet.ToArray();
                    int[] tempB = other.citySet.ToArray();
                    int size = tempA.Length;
                    for(int i = 0; i < size; i++) {
                        if (tempA[i] != tempB[i])
                            return false;
                    }
                    return true;
                }
                return false;
            }
        } 
        
        private Dictionary<MNode, double>[] memoize;
        private Dictionary<MNode, int>[] prev;
        #region Private members 

        /// <summary>
        /// Default number of cities (unused -- to set defaults, change the values in the GUI form)
        /// </summary>
        // (This is no longer used -- to set default values, edit the form directly.  Open Form1.cs,
        // click on the Problem Size text box, go to the Properties window (lower right corner), 
        // and change the "Text" value.)
        private const int DEFAULT_SIZE = 25;
        private const double INFINITY = double.PositiveInfinity;
        /// <summary>
        /// Default time limit (unused -- to set defaults, change the values in the GUI form)
        /// </summary>
        // (This is no longer used -- to set default values, edit the form directly.  Open Form1.cs,
        // click on the Time text box, go to the Properties window (lower right corner), 
        // and change the "Text" value.)
        private const int TIME_LIMIT = 60;        //in seconds

        private const int CITY_ICON_SIZE = 5;
        private Boolean initialBSSF = false;

        // For normal and hard modes:
        // hard mode only
        private const double FRACTION_OF_PATHS_TO_REMOVE = 0.20;

        /// <summary>
        /// the cities in the current problem.
        /// </summary>
        private City[] Cities;
        double[,] distance;
        /// <summary>
        /// a route through the current problem, useful as a temporary variable. 
        /// </summary>
        private ArrayList Route;
        /// <summary>
        /// best solution so far. 
        /// </summary>
        private TSPSolution bssf; 
        private int[] heldKarpPath;
        /// <summary>
        /// how to color various things. 
        /// </summary>
        private Brush cityBrushStartStyle;
        private Brush cityBrushStyle;
        private Pen routePenStyle;
        double currentBssf;

        private bool solutionFound;

        private PriorityQueue queue;
        private int maxStroedStates = 0;
        private int numBSSFUpdates = 0;
        private int totalStates = 0;
        private int numStatesPruned = 0;
        /// <summary>
        /// keep track of the seed value so that the same sequence of problems can be 
        /// regenerated next time the generator is run. 
        /// </summary>
        private int _seed;
        /// <summary>
        /// number of cities to include in a problem. 
        /// </summary>
        private int _size;

        /// <summary>
        /// Difficulty level
        /// </summary>
        private HardMode.Modes _mode;

        /// <summary>
        /// random number generator. 
        /// </summary>
        private Random rnd;

        /// <summary>
        /// time limit in milliseconds for state space search
        /// can be used by any solver method to truncate the search and return the BSSF
        /// </summary>
        private int time_limit;
        #endregion

        #region Public members

        /// <summary>
        /// These three constants are used for convenience/clarity in populating and accessing the results array that is passed back to the calling Form
        /// </summary>
        public const int COST = 0;           
        public const int TIME = 1;
        public const int COUNT = 2;
        
        public int Size
        {
            get { return _size; }
        }

        public int Seed
        {
            get { return _seed; }
        }
        #endregion

        #region Constructors
        public ProblemAndSolver()
        {
            this._seed = 1; 
            rnd = new Random(1);
            this._size = DEFAULT_SIZE;
            this.time_limit = TIME_LIMIT * 1000;                  // TIME_LIMIT is in seconds, but timer wants it in milliseconds

            this.resetData();
        }

        public ProblemAndSolver(int seed)
        {
            this._seed = seed;
            rnd = new Random(seed);
            this._size = DEFAULT_SIZE;
            this.time_limit = TIME_LIMIT * 1000;                  // TIME_LIMIT is in seconds, but timer wants it in milliseconds

            this.resetData();
        }

        public ProblemAndSolver(int seed, int size)
        {
            this._seed = seed;
            this._size = size;
            rnd = new Random(seed);
            this.time_limit = TIME_LIMIT * 1000;                        // TIME_LIMIT is in seconds, but timer wants it in milliseconds

            this.resetData();
        }
        public ProblemAndSolver(int seed, int size, int time)
        {
            this._seed = seed;
            this._size = size;
            rnd = new Random(seed);
            this.time_limit = time*1000;                        // time is entered in the GUI in seconds, but timer wants it in milliseconds

            this.resetData();
        }
        #endregion

        #region Private Methods

        /// <summary>
        /// Reset the problem instance.
        /// </summary>
        private void resetData()
        {

            Cities = new City[_size];
            Route = new ArrayList(_size);
            bssf = null;

            if (_mode == HardMode.Modes.Easy)
            {
                for (int i = 0; i < _size; i++)
                    Cities[i] = new City(rnd.NextDouble(), rnd.NextDouble());
            }
            else // Medium and hard
            {
                for (int i = 0; i < _size; i++)
                    Cities[i] = new City(rnd.NextDouble(), rnd.NextDouble(), rnd.NextDouble() * City.MAX_ELEVATION);
            }

            HardMode mm = new HardMode(this._mode, this.rnd, Cities);
            if (_mode == HardMode.Modes.Hard)
            {
                int edgesToRemove = (int)(_size * FRACTION_OF_PATHS_TO_REMOVE);
                mm.removePaths(edgesToRemove);
            }
            City.setModeManager(mm);

            cityBrushStyle = new SolidBrush(Color.Black);
            cityBrushStartStyle = new SolidBrush(Color.Red);
            routePenStyle = new Pen(Color.Blue,1);
            routePenStyle.DashStyle = System.Drawing.Drawing2D.DashStyle.Solid;
        }

        #endregion

        #region Public Methods

        /// <summary>
        /// make a new problem with the given size.
        /// </summary>
        /// <param name="size">number of cities</param>
        public void GenerateProblem(int size, HardMode.Modes mode)
        {
            this._size = size;
            this._mode = mode;
            resetData();
        }

        /// <summary>
        /// make a new problem with the given size, now including timelimit paremeter that was added to form.
        /// </summary>
        /// <param name="size">number of cities</param>
        public void GenerateProblem(int size, HardMode.Modes mode, int timelimit)
        {
            this._size = size;
            this._mode = mode;
            this.time_limit = timelimit*1000;                                   //convert seconds to milliseconds
            resetData();
        }

        /// <summary>
        /// return a copy of the cities in this problem. 
        /// </summary>
        /// <returns>array of cities</returns>
        public City[] GetCities()
        {
            City[] retCities = new City[Cities.Length];
            Array.Copy(Cities, retCities, Cities.Length);
            return retCities;
        }

        /// <summary>
        /// draw the cities in the problem.  if the bssf member is defined, then
        /// draw that too. 
        /// </summary>
        /// <param name="g">where to draw the stuff</param>
        public void Draw(Graphics g)
        {
            float width  = g.VisibleClipBounds.Width-45F;
            float height = g.VisibleClipBounds.Height-45F;
            Font labelFont = new Font("Arial", 10);

            // Draw lines
            if (bssf != null)
            {
                // make a list of points. 
                Point[] ps = new Point[bssf.Route.Count];
                int index = 0;
                foreach (City c in bssf.Route)
                {
                    if (index < bssf.Route.Count -1)
                        g.DrawString(" " + index +"("+c.costToGetTo(bssf.Route[index+1]as City)+")", labelFont, cityBrushStartStyle, new PointF((float)c.X * width + 3F, (float)c.Y * height));
                    else 
                        g.DrawString(" " + index +"("+c.costToGetTo(bssf.Route[0]as City)+")", labelFont, cityBrushStartStyle, new PointF((float)c.X * width + 3F, (float)c.Y * height));
                    ps[index++] = new Point((int)(c.X * width) + CITY_ICON_SIZE / 2, (int)(c.Y * height) + CITY_ICON_SIZE / 2);
                }

                if (ps.Length > 0)
                {
                    g.DrawLines(routePenStyle, ps);
                    g.FillEllipse(cityBrushStartStyle, (float)Cities[0].X * width - 1, (float)Cities[0].Y * height - 1, CITY_ICON_SIZE + 2, CITY_ICON_SIZE + 2);
                }

                // draw the last line. 
                g.DrawLine(routePenStyle, ps[0], ps[ps.Length - 1]);
            }

            // Draw city dots
            foreach (City c in Cities)
            {
                g.FillEllipse(cityBrushStyle, (float)c.X * width, (float)c.Y * height, CITY_ICON_SIZE, CITY_ICON_SIZE);
            }

        }

        /// <summary>
        ///  return the cost of the best solution so far. 
        /// </summary>
        /// <returns></returns>
        public double costOfBssf ()
        {
            if (bssf != null)
                return (bssf.costOfRoute());
            else
                return -1D; 
        }

        /// <summary>
        /// This is the entry point for the default solver
        /// which just finds a valid random tour 
        /// </summary>
        /// <returns>results array for GUI that contains three ints: cost of solution, time spent to find solution, number of solutions found during search (not counting initial BSSF estimate)</returns>
        public string[] defaultSolveProblem()
        {
            int i, swap, temp, count=0;
            string[] results = new string[3];
            int[] perm = new int[Cities.Length];
            Route = new ArrayList();
            Random rnd = new Random();
            Stopwatch timer = new Stopwatch();

            timer.Start();

            do
            {
                for (i = 0; i < perm.Length; i++)                                 // create a random permutation template
                    perm[i] = i;
                for (i = 0; i < perm.Length; i++)
                {
                    swap = i;
                    while (swap == i)
                        swap = rnd.Next(0, Cities.Length);
                    temp = perm[i];
                    perm[i] = perm[swap];
                    perm[swap] = temp;
                }
                Route.Clear();
                for (i = 0; i < Cities.Length; i++)                            // Now build the route using the random permutation 
                {
                    Route.Add(Cities[perm[i]]);
                }
                bssf = new TSPSolution(Route);
                count++;
            } while (costOfBssf() == double.PositiveInfinity);                // until a valid route is found
            timer.Stop();

            results[COST] = costOfBssf().ToString();                          // load results array
            results[TIME] = timer.Elapsed.ToString();
            results[COUNT] = count.ToString();

            return results;
        }

      
        /// <summary>
        /// performs a Branch and Bound search of the state space of partial tours
        /// stops when time limit expires and uses BSSF as solution
        /// </summary>
        /// <returns>results array for GUI that contains three ints: cost of solution, time spent to find solution, number of solutions found during search (not counting initial BSSF estimate)</returns>
        public string[] bBSolveProblem()
        {
            Stopwatch timer = new Stopwatch();
            timer.Start();
            initialBSSF = true;
            greedySolveProblem();            
            string[] results = new string[3];
            int n = Cities.Length;
            double[,] baseMatrix = new double[n,n];
            int count;
            ReducedMatrix currentCity;
            queue = new PriorityQueue();
            ReducedMatrix reduced;
            for(int i = 0; i < n; i++)
            {
                for(int j = 0; j < n; j++)
                {
                    if (i == j)
                        baseMatrix[i, j] = Double.PositiveInfinity;
                    else
                        baseMatrix[i, j] = Cities[i].costToGetTo(Cities[j]);
                }
            }
            reduced = new ReducedMatrix(baseMatrix, Cities.Length,Cities[0]);
            reduced.reduce();            
            queue.add(reduced);
            count = 0;
            totalStates = 1;
            while (true)
            {
                if(queue.getSize() > maxStroedStates)
                {
                    maxStroedStates = queue.getSize();
                }
                if (timer.ElapsedMilliseconds > time_limit)
                    break;
                currentCity = queue.remove();
                if (currentCity == null)
                    break;
                if (currentBssf < currentCity.bound)
                    break;
                for(int i=1;i< n; i++)
                {
                    if(currentCity.matrix[currentCity.cityNumber,i ] != INFINITY)
                    {
                        ReducedMatrix newMatrix = new ReducedMatrix(currentCity, Cities[i], i);
                        newMatrix.cancel(currentCity.cityNumber, i);
                        newMatrix.reduce();
                        if (newMatrix.Route.Count == Cities.Length)
                        {
                            if (newMatrix.bound < currentBssf)
                            {
                                numBSSFUpdates++;
                                bssf = new TSPSolution(newMatrix.Route);
                                currentBssf = newMatrix.bound;
                            }
                            count++;
                            break;
                        }
                        if (newMatrix.bound < currentBssf)
                        {                            
                            queue.add(newMatrix);
                            totalStates++;                     
                        }
                        else
                        {
                            numStatesPruned++;
                        }
                    }
                }
            }

            timer.Stop();
            numStatesPruned += queue.getSize();

            String track;
            track = "MAX #of stored states " + maxStroedStates + "\n";
            track += "#of BSSF updates " + numBSSFUpdates +"\n";
            track += "Total states created " + totalStates + "\n";
            track += "Total states pruned " + numStatesPruned + "\n";

            MessageBox.Show(track);


            results[COST] = bssf.costOfRoute().ToString();    // load results into array here, replacing these dummy values
            results[TIME] = timer.Elapsed.ToString();
            results[COUNT] = count.ToString() ;

            return results;
        }

        /////////////////////////////////////////////////////////////////////////////////////////////
        // These additional solver methods will be implemented as part of the group project.
        ////////////////////////////////////////////////////////////////////////////////////////////

        /// <summary>
        /// finds the greedy tour starting from each city and keeps the best (valid) one
        /// </summary>
        /// <returns>results array for GUI that contains three ints: cost of solution, time spent to find solution, number of solutions found during search (not counting initial BSSF estimate)</returns>
        public string[] greedySolveProblem()
        {
            Stopwatch timer = new Stopwatch();

            timer.Start();
            string[] results = new string[3];
            Boolean[] visited = new Boolean[Cities.Length];
            double min;
            int temp;
            int currentCity;
            double cost;
            double totalCost;
            double minCost = INFINITY;
            int count=0;
            for (int i = 0; i < Cities.Length; i++)
            {
                if (timer.Elapsed.TotalMilliseconds > time_limit)
                    break;
                for(int j = 0; j < Cities.Length; j++)
                {
                    visited[j] = false;
                }
                visited[i] = true;
                Route = new ArrayList();
                Route.Add(Cities[i]);
                currentCity = i;
                totalCost = 0;
                temp = -1;
                for (int j=0; j < Cities.Length - 1; j++)
                {
                    min = INFINITY;                    
                    for(int k=0; k< Cities.Length; k++)
                    {
                        if (visited[k] || i==k)
                            continue;
                        cost = Cities[currentCity].costToGetTo(Cities[k]);
                        if (cost == double.PositiveInfinity)
                            continue;
                        if (min > Cities[currentCity].costToGetTo(Cities[k]))
                        {
                            min = Cities[currentCity].costToGetTo(Cities[k]);
                            temp = k;
                        }
                    }
                    totalCost += min;
                    visited[temp] = true;
                    Route.Add(Cities[temp]);
                }
                count++;
                totalCost += Cities[temp].costToGetTo(Cities[i]);               
                if( totalCost < minCost && totalCost != INFINITY)
                {
                    bssf = new TSPSolution(Route);
                    minCost = totalCost;
                    if (initialBSSF)
                    {
                        currentBssf = minCost;
                        return null;
                    }

                    results[COST] = minCost.ToString();    // load results into array here, replacing these dummy values
                    results[TIME] = timer.Elapsed.ToString();
                    results[COUNT] = count.ToString();

                    return results;
                }
            }
            timer.Stop();

            results[COST] = minCost.ToString();    // load results into array here, replacing these dummy values
            results[TIME] = timer.Elapsed.ToString();
            results[COUNT] = count.ToString();

            return results;
        }
        
        private double getCost(HashSet<int> set, int prevVertex, int setSize) {
            HashSet<int> tempSet = new HashSet<int>(set);
            tempSet.Remove(prevVertex);
            double cost = double.MaxValue;
            MNode index = new MNode { Vertex = prevVertex, citySet = tempSet };
            if (setSize == 0) {
                return distance[0, prevVertex];
            }
            else {
                Dictionary<MNode, double> dictionary = memoize[tempSet.Count];
                cost = dictionary[index];
            }
            /*
            for(int i = 0; i < memoize[setSize-1].Count; i++) {
                if(memoize[setSize-1].ElementAt(i).Key.Equals(index)) {
                    cost = memoize[setSize - 1].ElementAt(i).Value;
                    break;
                }
            }    */        
            return cost;
        }

        public List<HashSet<int>> gemerateCombination(int n) {
            int[] input = new int[n];
            List<HashSet<int>> allSets = new List<HashSet<int>>();
            for (int i = 0; i < n; i++) {
                generateCombinations(n, 0, 1, i, new HashSet<int>(), ref allSets);
            }            
            return allSets;
        }

        public void generateCombinations(int n, int pos, int index, int size, HashSet<int> currentSet, ref List<HashSet<int>> allSets) {
            if(pos == size) {                
                allSets.Add(currentSet);
                return;
            }
            if (index == n)
                return;
            // insert index into the set
            currentSet.Add(index);
            generateCombinations(n, pos + 1, index + 1, size, new HashSet<int>(currentSet), ref allSets);
            // not insertig index into the set
            currentSet.Remove(index);
            generateCombinations(n, pos, index + 1, size, currentSet, ref allSets);
        }

        private double GetMinimumCostRoute(int n) {
            List<HashSet<int>> allSets = gemerateCombination(n);
            for(int i = 0; i < n; i++) {
                memoize[i] = new Dictionary<MNode, double>();
                prev[i] = new Dictionary<MNode, int>();
            }
            foreach (HashSet<int> set in allSets) {
                for (int currentVertex = 1; currentVertex < n; currentVertex++) {
                    if (set.Contains(currentVertex))
                        continue;
                    MNode index = new MNode { Vertex = currentVertex, citySet = set };                    
                    int minPrevVertex = -1;
                    //to avoid ConcurrentModificationException copy set into another set while iterating
                    HashSet<int> copySet = new HashSet<int>(set);
                    double minCost = double.MaxValue;
                    foreach (int prevVertex in set) {
                        if(distance[prevVertex, currentVertex] == INFINITY) {
                            continue;
                        }
                        double cost = distance[prevVertex,currentVertex] + getCost(copySet, prevVertex, set.Count);
                        if (cost < minCost) {
                            minCost = cost;
                            minPrevVertex = prevVertex;
                        }
                    }
                    //this happens for empty subset
                    if (set.Count == 0) {
                        minCost = distance[0,currentVertex];
                    }
                    memoize[set.Count].Add(index, minCost);
                    prev[set.Count].Add(index, minPrevVertex);
                }
            }

            HashSet<int> tempset = new HashSet<int>();
            for (int i = 1; i < n; i++) {
                tempset.Add(i);
            }
            double min = INFINITY;
            int prevvertex = -1;
            //to avoid ConcurrentModificationException copy set into another set while iterating
            HashSet<int> copyset = new HashSet<int>(tempset);
            double endCost;
            foreach (int l in tempset) {
                endCost = distance[l,0] + getCost(copyset, l, n-1);
                if (endCost < min) {
                    min = endCost;
                    prevvertex = l;
                }
            }

            MNode finalNode = new MNode { Vertex = 0, citySet = tempset };
            memoize[n-1].Add(finalNode, min);
            double totalCost = min;
            prev[n-1].Add(finalNode, prevvertex);

            ArrayList tour = new ArrayList();            
            int currentCity = 0;
            int k = n - 1;
            for (int i = 0; i < n; i++) {
                tour.Add(Cities[currentCity]);
                MNode node = new MNode { Vertex = currentCity, citySet = tempset };
                currentCity = prev[k][node];
                /*for (int j = 0; j < memoize[k].Count; j++) {
                    if (prev[k].ElementAt(j).Key.Equals(node)) {
                        currentCity = prev[k].ElementAt(j).Value;
                        break;
                    }
                }*/
                k--;
                tempset.Remove(currentCity);
            }
            tour.Reverse();

            bssf = new TSPSolution(tour);
            return totalCost;
        }


        // implementation of Held-Karp Algorithm
        public string[] fancySolveProblem()
        {
            string[] results = new string[3];            
            Stopwatch timer = new Stopwatch();
            timer.Start();
            initialBSSF = true;
            greedySolveProblem();
            solutionFound = false;
            int n = Cities.Length;
            distance = new Double[n, n];
            var set = new HashSet<int>();            
            heldKarpPath = new int[n];
            memoize = new Dictionary<MNode, double>[n];
            prev = new Dictionary<MNode, int>[n];
            currentBssf = double.MaxValue;
            for (int i = 0; i < n; i++) {
                set.Add(i);
                for (int j = 0; j < n; j++) {
                    if (i == j)
                        distance[i, j] = Double.PositiveInfinity;
                    else
                        distance[i, j] = Cities[i].costToGetTo(Cities[j]);
                }
            }
           set.Remove(0);
            double cost = GetMinimumCostRoute(n);
            MNode node = new MNode { Vertex = 0, citySet = set };
            
            timer.Stop();
            
            results[COST] = bssf.costOfRoute().ToString();    // load results into array here, replacing these dummy values
            results[TIME] = timer.Elapsed.ToString();
            results[COUNT] = "1";

            return results;
        }
        #endregion
    }
}
