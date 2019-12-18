import java.util.Collections;
import java.util.List;
import java.util.ArrayList;

public class TspDynamic {

    private int N, start;
    private final double[][] cost;
    private List<Integer> path = new ArrayList<>();
    private double minPathCost = Double.MAX_VALUE;
    private boolean ranSolver = false;

    // constructor to get cost matrix
    public TspDynamic (double[][] cost) {
        this(0, cost);
    }

    // constructor including starting node. Can make starting node any node.
    public TspDynamic (int start, double[][] cost) {
        N = cost.length;
        if (start < 0 || start >= N) throw new IllegalArgumentException("Invalid start node.");

        this.start = start;
        this.cost = cost;
    }

    // Method that returns the optimal tour for the traveling salesman problem
    public List<Integer> getPath() {
        if (!ranSolver)
            solve();
        return path;  // returns the optimal path
    }

    // Return the minimal path cost
    public double getPathCost() {
        if (!ranSolver)
            solve();
        return minPathCost; // return the minimum cost
    }

    // Solves the traveling salesman problem and catches solution
    public void solve() {
        if (ranSolver)
            return;
        final int END_STATE = (1 << N) - 1;  // Where all node have been visited and set it to 1

        // Initialize table to null
        Double[][] memo = new Double[N][1 << N];

        // Add all outgoing edges from the starting node to memo table
        for (int end = 0; end < N; ++end) {
            if (end == start)
                continue;
            memo[end][(1 << start) | (1 << end)] = cost[start][end];
        }

        for (int node = 3; node <= N; ++node) {
            for (int subset : combinations(node, N)) {
                // Need to make sure it exists in the subset
                if (notIn(start, subset))
                    continue;
                // Here we need to make sure the next node isn't the starting node and it's in subset
                for (int next = 0; next < N; ++next) {
                    if (next == start || notIn(next, subset))
                        continue;
                    int subsetWithoutNext = subset ^ (1 << next);  // make next node 0
                    double minCost = Double.MAX_VALUE;
                    // For every end node that isn't the start or end and exists in subset
                    // calculate the new cost
                    for (int end = 0; end < N; ++end) {
                        // if the next cost is less than minCost then record that in memo
                        if (end == start || end == next || notIn(end, subset))
                            continue;
                        double newCost = memo[end][subsetWithoutNext] + cost[end][next];
                        if (newCost < minCost) {
                            minCost = newCost;
                        }
                    }
                    memo[next][subset] = minCost;
                }
            }
        }

        // Connect path back to starting node and minimize cost
        for (int i = 0; i < N; ++i) {
            if (i == start)
                continue;
            double pathCost = memo[i][END_STATE] + cost[i][start];
            if (pathCost < minPathCost) {
                minPathCost = pathCost;  // update the minPathCost which shrinks it from MAX_VALUE
            }
        }

        // FIND THE PATH
        int lastIndex = start;
        int state = END_STATE;
        path.add(start);

        // Reconstruct TSP path from memo table
        for (int i = 1; i < N; ++i) {
            int index = -1;
            // Find the next best node
            for (int j = 0; j < N; ++j) {
                if (j == start || notIn(j, state))
                    continue;  // Skip starting node
                if (index == -1)
                    index = j;
                // Look in table for end node i and the END_STATE so the path goes back to start
                // Check which is smaller
                double preCost = memo[index][state] + cost[index][lastIndex];
                double newCost = memo[j][state] + cost[j][lastIndex];
                // Check if the newCost is less than the previous one
                if (newCost < preCost) {
                    index = j;
                }
            }
            // Add it to our path and set state to zero and make lastIndex to index
            path.add(index);
            state = state ^ (1 << index);
            lastIndex = index;
        }

        // Add that starting node to path and reverse the order
        path.add(start);
        Collections.reverse(path);

        // Set our boolean to true
        ranSolver = true;
    }

    // Check if the element is set
    private static boolean notIn (int elem, int subset) {
        return ((1 << elem) & subset) == 0;
    }

    // Method to set nodes to 1 if bit sets are found of size n
    public static List<Integer> combinations (int node, int n) {
        List<Integer> subsets = new ArrayList<>();
        combinations(0, 0, node, n, subsets);
        return subsets;
    }

    // Find all the combinations of size node.
    private static void combinations(int set, int at, int node, int n, List<Integer> subsets) {
        // Return if there are more elements left to select than what is available
        int elementsLeftToPick = n - at;
        if (elementsLeftToPick < node)
            return;

        // Execute if node has the subset
        if (node == 0) {
            subsets.add(set);
        }
        else {
            for (int i = at; i < n; ++i) {
                // Try to add this
                set |= 1 << i;
                combinations(set, i + 1, node - 1, n, subsets);

                // Go back and try the instance where we did not include this element
                // BTW I had to learn some bitwise operators for this
                // https://www.tutorialspoint.com/cprogramming/c_operators.htm
                set &= ~(1 << i);
            }
        }
    }

    public static void main (String[] args) {
        int n = (int)(10 * Math.random());

        // make sure n > 2
        while (n <= 2)
            n = (int)(10 * Math.random());

        System.out.println("My n: " + n);
        double[][] costMatrix = new double[n][n];

        // Create adjacency matrix
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                costMatrix[i][j] = (int)(10 * Math.random());
            }
        }



        /*
        // This matrix for correctness testing

        int n = 5;
        double[][] costMatrix = new double[n][n];
        for (double[] row : costMatrix) java.util.Arrays.fill(row, 100);
        costMatrix[1][0] = 1;
        costMatrix[2][1] = 1;
        costMatrix[3][2] = 1;
        costMatrix[4][3] = 1;
        costMatrix[0][4] = 1;
         */


        int startNode = 0;
        TspDynamic solver = new TspDynamic(startNode, costMatrix);
        System.out.println("Path: " + solver.getPath());
        System.out.println("Cost: " + solver.getPathCost());
    }
}
