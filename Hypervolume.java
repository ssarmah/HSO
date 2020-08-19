import java.io.IOException;
import java.util.concurrent.TimeUnit;
/**
 *
 * @author IIITG
 */
public class Hypervolume {
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException {
        // TODO code application logic here
        Helper helper = new Helper();
        
		int noSol = 5; // number of rows (solutions)
		int noCol = 4; // number of columns (objectives)
		/* Initialize the population */
		Global.population = new double[noSol][noCol];
		Global.population[0][0] = 5; Global.population[0][1] = 5; Global.population[0][2] = 5; Global.population[0][3] = 1;
		Global.population[1][0] = 4; Global.population[1][1] = 4; Global.population[1][2] = 4; Global.population[1][3] = 2;
		Global.population[2][0] = 3; Global.population[2][1] = 3; Global.population[2][2] = 3; Global.population[2][3] = 3;
		Global.population[3][0] = 2; Global.population[3][1] = 2; Global.population[3][2] = 2; Global.population[3][3] = 4;
		Global.population[4][0] = 1; Global.population[4][1] = 1; Global.population[4][2] = 1; Global.population[4][3] = 5;
        
		/* Fix the reference point */
        Global.refPoint = new double[noCol];
        for(int m = 0; m < noCol; m++) {
            Global.refPoint[m] = 0;
        }
			
        /* Assign each solution an id; i-th solution has id i */
        int populationId[] = new int[noSol];
        for (int i = 0; i < noSol; i++) { 
            populationId[i] = i;
        }

        /* Initially all the objectives are considered starting from first to last */
        int lb = 0; // The first objective is '0'-th
        int ub = noCol-1; // The last objective is 'noCol-1'-th

		/* Set number of explored node to be 0 */
        Global.noNodes_HSO = 0; // By HSO
        Global.noNodes_IHSO = 0; // By improved HSO

        System.out.println("The initial population is ...");
        helper.print(populationId, lb, ub); // print the solutions considering all the objectives

        /* Calculate the hypervolume by HSO */
        long startTime = System.nanoTime();
        double volume1 = helper.calculateHyperVolume_Version1(populationId, lb, ub);
        System.out.println("Hypervolume = " + volume1);
        System.out.println("Number of nodes exploed by HSO = " + (Global.noNodes_HSO+1));
        long endTime = System.nanoTime();
        long durationInNano_HSO = endTime - startTime; 
        System.out.println("Time in nano seconds by HSO: " + durationInNano_HSO);
                
                
                
                
        /* Calculate the hypervolume by IHSO */
       startTime = System.nanoTime();
        Global.hashMatrix = new Trie[noCol][noSol];
        for(int i = 0; i < noCol; i++) {
            for(int j = 0; j < noSol; j++) {
                Global.hashMatrix[i][j] = new Trie();
            }
        }
        double volume3 = helper.calculateHyperVolume_Version3(populationId, lb, ub);
        System.out.println("Hypervolume = " + volume3);
        System.out.println("Number of nodes exploed by IHSO = " + (Global.noNodes_IHSO+1));
        endTime = System.nanoTime();
        long durationInNano_IHSO = endTime - startTime; 
        System.out.println("Time in nano seconds by IHSO: " + durationInNano_IHSO);
    } 
}



class Global {
    public static double[][] population;  
    public static double[] refPoint;
    
    public static double noNodes_HSO = 0;
    public static double noNodes_IHSO = 0;
    
    public static Trie[][] hashMatrix;    
}

class HeapSort {
    /* Code is inspired from https://www.geeksforgeeks.org/heap-sort-for-decreasing-order-using-min-heap/ */
    // To heapify a subtree rooted with node i which is
    // an index in arr[]. n is size of heap
    // To heapify a subtree rooted with node i which is 
    // an index in arr[]. n is size of heap 
    void heapify(int arr[], int n, int i, int obj) { 
        int smallest = i; // Initialize smallest as root 
        int l = 2 * i + 1; // left = 2*i + 1 
        int r = 2 * i + 2; // right = 2*i + 2 
  
        // If left child is smaller than root 
        if (l < n && Global.population[arr[l]][obj] < Global.population[arr[smallest]][obj]) 
            smallest = l; 
  
        // If right child is smaller than smallest so far 
        if (r < n && Global.population[arr[r]][obj] < Global.population[arr[smallest]][obj]) 
            smallest = r; 
  
        // If smallest is not root 
        if (smallest != i) { 
            int temp = arr[i]; 
            arr[i] = arr[smallest]; 
            arr[smallest] = temp; 
  
            // Recursively heapify the affected sub-tree 
            heapify(arr, n, smallest, obj); 
        } 
    } 
    
    /* Function to sort the solutions based on obj-th objective */
    int[] heapSort(int arr[], int obj) { 
        // Build heap (rearrange array)
        int n = arr.length;
        for (int i = n / 2 - 1; i >= 0; i--) 
            heapify(arr, n, i, obj); 
  
        // One by one extract an element from heap 
        for (int i = n - 1; i >= 0; i--) { 
              
            // Move current root to end 
            int temp = arr[0]; 
            arr[0] = arr[i]; 
            arr[i] = temp; 
  
            // call max heapify on the reduced heap 
            heapify(arr, i, 0, obj); 
        } 
        return arr;
    } 
}


class Helper {
    
    /* Function to compute the hypervolume of a set of solutions in a recursive manner.
        Base case: Number of objectives = 1 */
    public double calculateHyperVolume_Version1(int[] populationId, int lb, int ub) {
        /*System.out.println("\nThe population inside calculateHyperVolume is ...");
        print(populationId, lb, ub);*/
        
        int[] nonDominatedPopulationId;
        nonDominatedPopulationId = obtainNonDominatedPopulation(populationId, lb, ub);
        
        /*System.out.println("\nThe non-dominated population inside calculateHyperVolume is ...");
        print(nonDominatedPopulationId, lb, ub);*/
        
        HeapSort hs = new HeapSort();
        /* Sort the non-dominated population based on lb-th objective */
        nonDominatedPopulationId = hs.heapSort(nonDominatedPopulationId, lb);
        /*System.out.println("\nThe Sorted non-dominated population after sorting inside calculateHyperVolume is ...");
        print(nonDominatedPopulationId, lb, ub);*/
        
        if (lb == ub) { // Only one objective is left
            /* In case one objective left and we have one solution, 
                then the value of the solution in that particular will be the volume
                
               In case one objective left and we have more than one solution, 
                then the value of the first solution in that particular will be the volume. 
                We have alrday sorted the solutions so the first solution will dominate others in case 
                all the solutions have unique value for lb-th objective. Even if the solutions have the same value for the 
                lb-th objective, then any value can be considered (as all are same and we have considered the first)
            */
            //System.out.println("     volume = " + Global.population[nonDominatedPopulationId[0]][lb]);
            return Global.population[nonDominatedPopulationId[0]][lb];
        }
        
        
        double volume = 0;
        for (int i = 0; i < nonDominatedPopulationId.length; i++) { // For each slice
            Global.noNodes_HSO++;
            /* Obtain the slice */
            int[] slice = new int[i+1];
            for (int j = 0; j < i+1; j++) {
                slice[j] = nonDominatedPopulationId[j];
            }
            
            /* Recursively compute the volume of the slice considering next objective (lb+1) to create the slice */
            double sliceVol = calculateHyperVolume_Version1(slice, lb+1, ub); 
            //System.out.println("sliceVol = " + sliceVol);
            /* Find the depth to which the volume of the slice gets multiplied */
            double depth;
            if(i == nonDominatedPopulationId.length-1) {
                depth = Global.population[nonDominatedPopulationId[i]][lb] - Global.refPoint[lb];
            } else {
                depth = Global.population[nonDominatedPopulationId[i]][lb] - Global.population[nonDominatedPopulationId[i+1]][lb];
            }  
            volume = volume + sliceVol * depth;
        }
        //System.out.println("volume = " + volume);
        return volume;
    }
    
    
    
    /* Function to compute the hypervolume of a set of solutions in a recursive manner 
        Base case: Number of solutions/objectives = 1 */
    public double calculateHyperVolume_Version2(int[] populationId, int lb, int ub) {
        /*System.out.println("\nThe population inside calculateHyperVolume is ...");
        print(populationId, lb, ub);*/
        
        if(populationId.length == 1) {
            double volume = 1;
            for (int j = lb; j <= ub; j++) {
                volume = volume * Global.population[populationId[0]][j]; 
            }
            return volume;
        }
        
        int[] nonDominatedPopulationId;
        nonDominatedPopulationId = obtainNonDominatedPopulation(populationId, lb, ub);
        
        /*System.out.println("\nThe non-dominated population inside calculateHyperVolume is ...");
        print(nonDominatedPopulationId, lb, ub);*/
        
        HeapSort hs = new HeapSort();
        /* Sort the non-dominated population based on lb-th objective */
        nonDominatedPopulationId = hs.heapSort(nonDominatedPopulationId, lb);
        /*System.out.println("\nThe Sorted non-dominated population after sorting inside calculateHyperVolume is ...");
        print(nonDominatedPopulationId, lb, ub);*/
        
        if (lb == ub) { // Only one objective is left
            /* In case one objective left and we have one solution, 
                then the value of the solution in that particular will be the volume
                
               In case one objective left and we have more than one solution, 
                then the value of the first solution in that particular will be the volume. 
                We have alrday sorted the solutions so the first solution will dominate others in case 
                all the solutions have unique value for lb-th objective. Even if the solutions have the same value for the 
                lb-th objective, then any value can be considered (as all are same and we have considered the first)
            */
            return Global.population[nonDominatedPopulationId[0]][lb];
        }

        double volume = 0;
        for (int i = 0; i < nonDominatedPopulationId.length; i++) { // For each slice
            /* Obtain the slice */
            int[] slice = new int[i+1];
            for (int j = 0; j < i+1; j++) {
                slice[j] = nonDominatedPopulationId[j];
            }
            
            /* Recursively compute the volume of the slice considering next objective (lb+1) to create the slice */
            double sliceVol = calculateHyperVolume_Version2(slice, lb+1, ub); 
            
            /* Find the depth to which the volume of the slice gets multiplied */
            double depth;
            if(i == nonDominatedPopulationId.length-1) {
                depth = Global.population[nonDominatedPopulationId[i]][lb] - Global.refPoint[lb];
            } else {
                depth = Global.population[nonDominatedPopulationId[i]][lb] - Global.population[nonDominatedPopulationId[i+1]][lb];
            }  
            
            volume = volume + sliceVol * depth;
        }
        return volume;
    }
    
    
    
    /* Function to compute the hypervolume of a set of solutions in a recursive manner 
        Base case: Number of solutions/objectives = 1 */
    public double calculateHyperVolume_Version3(int[] populationId, int lb, int ub) {
        /*System.out.println("\nThe population inside calculateHyperVolume is ...");
        print(populationId, lb, ub);*/
        
        int aux[] = countingSort(populationId);
        //System.out.println("Aux array: " + Arrays.toString(aux));
        
        double value = Global.hashMatrix[ub-lb][populationId.length-1].search(aux);
        //System.out.println("value = " + value);
        
        
        if(value == -1) { // populationId is not alreday processed
            //System.out.println("Not Found");
            /*if(populationId.length == 1) {
                double volume = 1;
                for (int j = lb; j <= ub; j++) {
                    volume = volume * Global.population[populationId[0]][j]; 
                }
                
                Global.hashMatrix[ub-lb][populationId.length-1].insert(aux, volume);
                print();
                return volume;
            }*/

            int[] nonDominatedPopulationId;
            nonDominatedPopulationId = obtainNonDominatedPopulation(populationId, lb, ub);

            /*System.out.println("\nThe non-dominated population inside calculateHyperVolume is ...");
            print(nonDominatedPopulationId, lb, ub);*/

            HeapSort hs = new HeapSort();
            /* Sort the non-dominated population based on lb-th objective */
            nonDominatedPopulationId = hs.heapSort(nonDominatedPopulationId, lb);
            /*System.out.println("\nThe Sorted non-dominated population after sorting inside calculateHyperVolume is ...");
            print(nonDominatedPopulationId, lb, ub);*/

            if (lb == ub) { // Only one objective is left
                /* In case one objective left and we have one solution, 
                    then the value of the solution in that particular will be the volume

                   In case one objective left and we have more than one solution, 
                    then the value of the first solution in that particular will be the volume. 
                    We have alrday sorted the solutions so the first solution will dominate others in case 
                    all the solutions have unique value for lb-th objective. Even if the solutions have the same value for the 
                    lb-th objective, then any value can be considered (as all are same and we have considered the first)
                */
                double volume = Global.population[nonDominatedPopulationId[0]][lb];
                //System.out.println("volume (m=1) = " + volume);
                Global.hashMatrix[ub-lb][populationId.length-1].insert(aux, volume);
                return volume;
            }


            double volume = 0;
            for (int i = 0; i < nonDominatedPopulationId.length; i++) { // For each slice
                Global.noNodes_IHSO++;
                /* Obtain the slice */
                int[] slice = new int[i+1];
                for (int j = 0; j < i+1; j++) {
                    slice[j] = nonDominatedPopulationId[j];
                }

                /* Recursively compute the volume of the slice considering next objective (lb+1) to create the slice */
                double sliceVol = calculateHyperVolume_Version3(slice, lb+1, ub); 

                /* Find the depth to which the volume of the slice gets multiplied */
                double depth;
                if(i == nonDominatedPopulationId.length-1) {
                    depth = Global.population[nonDominatedPopulationId[i]][lb] - Global.refPoint[lb];
                } else {
                    depth = Global.population[nonDominatedPopulationId[i]][lb] - Global.population[nonDominatedPopulationId[i+1]][lb];
                }  

                volume = volume + sliceVol * depth;
            }
            //System.out.println("volume = " + volume);
            Global.hashMatrix[ub-lb][populationId.length-1].insert(aux, volume);
            //print();
            return volume;
        } else {
            //System.out.println("****************************");
            return value;
        }
    }
    
    
    
    
    
    
    
    /* Function to compute the hypervolume of a set of solutions in a recursive manner 
        Base case: Number of solutions/objectives = 1 */
    public double calculateHyperVolume_Version4(int[] populationId, int lb, int ub) {
        /*System.out.println("\nThe population inside calculateHyperVolume is ...");
        print(populationId, lb, ub);*/
        
        
        if(populationId.length == 1) {
            double volume = 1;
            for (int j = lb; j <= ub; j++) {
                volume = volume * Global.population[populationId[0]][j]; 
            }
            //System.out.println("volume (n=1) = " + volume);
            // Do not add in HASH when we have single solution in slice
            //print();
            return volume;
       }
        
        
        int aux[] = countingSort(populationId);
        //System.out.println("Aux array: " + Arrays.toString(aux));
        
        double value = Global.hashMatrix[ub-lb][populationId.length-1].search(aux);
        //System.out.println("value = " + value);
        
        
        if(value == -1) { // populationId is not alreday processed
            //System.out.println("Not Found");
            /*if(populationId.length == 1) {
                double volume = 1;
                for (int j = lb; j <= ub; j++) {
                    volume = volume * Global.population[populationId[0]][j]; 
                }
                System.out.println("volume (n=1) = " + volume);
                
                Global.hashMatrix[ub-lb][populationId.length-1].insert(aux, volume);
                print();
                return volume;
            }*/

            int[] nonDominatedPopulationId;
            nonDominatedPopulationId = obtainNonDominatedPopulation(populationId, lb, ub);

            /*System.out.println("\nThe non-dominated population inside calculateHyperVolume is ...");
            print(nonDominatedPopulationId, lb, ub);*/

            HeapSort hs = new HeapSort();
            /* Sort the non-dominated population based on lb-th objective */
            nonDominatedPopulationId = hs.heapSort(nonDominatedPopulationId, lb);
            /*System.out.println("\nThe Sorted non-dominated population after sorting inside calculateHyperVolume is ...");
            print(nonDominatedPopulationId, lb, ub);*/

            if (lb == ub) { // Only one objective is left
                /* In case one objective left and we have one solution, 
                    then the value of the solution in that particular will be the volume

                   In case one objective left and we have more than one solution, 
                    then the value of the first solution in that particular will be the volume. 
                    We have alrday sorted the solutions so the first solution will dominate others in case 
                    all the solutions have unique value for lb-th objective. Even if the solutions have the same value for the 
                    lb-th objective, then any value can be considered (as all are same and we have considered the first)
                */
                double volume = Global.population[nonDominatedPopulationId[0]][lb];
                //System.out.println("volume (m=1) = " + volume);
                Global.hashMatrix[ub-lb][populationId.length-1].insert(aux, volume);
                return volume;
            }


            double volume = 0;
            for (int i = 0; i < nonDominatedPopulationId.length; i++) { // For each slice
                /* Obtain the slice */
                int[] slice = new int[i+1];
                for (int j = 0; j < i+1; j++) {
                    slice[j] = nonDominatedPopulationId[j];
                }

                /* Recursively compute the volume of the slice considering next objective (lb+1) to create the slice */
                double sliceVol = calculateHyperVolume_Version4(slice, lb+1, ub); 

                /* Find the depth to which the volume of the slice gets multiplied */
                double depth;
                if(i == nonDominatedPopulationId.length-1) {
                    depth = Global.population[nonDominatedPopulationId[i]][lb] - Global.refPoint[lb];
                } else {
                    depth = Global.population[nonDominatedPopulationId[i]][lb] - Global.population[nonDominatedPopulationId[i+1]][lb];
                }  

                volume = volume + sliceVol * depth;
            }
            //System.out.println("volume = " + volume);
            Global.hashMatrix[ub-lb][populationId.length-1].insert(aux, volume);
            //print();
            return volume;
        } else {
            //System.out.println("****************************");
            return value;
        }
    }
    
    
    
    public int[] obtainNonDominatedPopulation(int[] populationId, int lb, int ub) {
        int noSol = populationId.length;
        /* Number of non-dominated solutions in the population */
        int noNonDominatedSol = 0;
        
        boolean[] isDominated = new boolean[noSol];
        
        for(int i = 0; i < noSol; i++) {
            if(isDominated[i] == false) {
                boolean isSol_i_Dominated = false; // Solution i is not dominated
                for(int j = i+1; j < noSol; j++) {
                    int domRelation = dominanceRelation(populationId[i], populationId[j], lb, ub);
                    //System.out.println("domRelation = " + domRelation);
                    if(domRelation == 1) { // population[i] dominates population[j]
                        isDominated[j] = true; 
                    } else if (domRelation == -1) { // population[j] dominates population[i]
                        isDominated[i] = true;
                        isSol_i_Dominated = true;
                        break;
                    } 
                }
                if(isSol_i_Dominated == false) { // Solution i is not dominated
                    noNonDominatedSol++;
                }
            }
        }

        /* Collect all non-dominated solutions */
        int[] nonDominatedPopulationId = new int[noNonDominatedSol];
        int index = 0;
        for(int i = 0; i < noSol; i++) {
            if(isDominated[i] == false) {
                nonDominatedPopulationId[index] = populationId[i];
                index++;
            }
        }
        return nonDominatedPopulationId;
    }
    
    /* Obtain the dominance relationship between two solutions: sol1 and sol2
        considering lb-th to ub-th objective */
    public int dominanceRelation(int sol1, int sol2, int lb, int ub){
        boolean sol1_dom_sol2 = false;
        boolean sol2_dom_sol1 = false;
        for(int i = lb; i <= ub; i++) {
            if(Global.population[sol1][i] > Global.population[sol2][i]) {
                if(sol2_dom_sol1 == true) 
                    return 0;
                if(sol1_dom_sol2 == false) 
                    sol1_dom_sol2 = true;     
            } else if(Global.population[sol1][i] < Global.population[sol2][i]) {
                if(sol1_dom_sol2 == true) 
                    return 0;
                if(sol2_dom_sol1 == false) 
                    sol2_dom_sol1 = true;    
            }
        }
        if(sol1_dom_sol2 == true) {
            return 1;
        } else if(sol2_dom_sol1 == true) {
            return -1;
        } else {
            return 0;
        }
    }
    
    
    public void print(int[] populationId, int lb, int ub) {
        for (int i = 0; i < populationId.length; i++) { 
            System.out.print("Solution-" + i + " = [");
            for(int j = lb; j <= ub; j++) {
                if (j == ub)
                    System.out.print(Global.population[populationId[i]][j]); 
                else 
                    System.out.print(Global.population[populationId[i]][j] + ", "); 
            }
            System.out.print("]");
            System.out.println();
        }
    }
    
    
    
    public void print() {
        int noCol = Global.hashMatrix.length;
        int noSol = Global.hashMatrix[0].length;
        
        System.out.println("Global.hashMatrix ..."); 
        for (int i = 0; i < noCol; i++) { 
            for (int j = 0; j < noSol; j++) { 
                if (Global.hashMatrix[i][j].root == null) {
                    System.out.print("NULL  "); 
                } else {
                    System.out.print("****  "); 
                }
            }
            System.out.println(); 
        }
        System.out.println();        
    }
    
    
    public int[] countingSort(int[] array) { 
        int[] aux = new int[array.length];

        // find the smallest and the largest value
        int min = array[0];
        int max = array[0];
        for (int i = 1; i < array.length; i++) {
            if (array[i] < min) 
                min = array[i];
            else if (array[i] > max) 
                max = array[i];
        }

        // init array of frequencies
        int[] counts = new int[max - min + 1];

        // init the frequencies
        for (int i = 0;  i < array.length; i++)
          counts[array[i] - min]++;

        // recalculate the array - create the array of occurences
        counts[0]--;
        for (int i = 1; i < counts.length; i++) 
          counts[i] = counts[i] + counts[i-1];

        for (int i = array.length - 1; i >= 0; i--) 
            aux[counts[array[i] - min]--] = array[i];
        
        int[] auxDescending = new int[array.length];
        for (int i = 0; i < array.length; i++) {
            auxDescending[i] = aux[array.length - i - 1];
        }
        //return auxDescending;
        return aux;
    } 
}


class Trie {
    // Root of BST 
    TrieNode root; 
  
    // Constructor 
    Trie() {  
        this.root = null; 
    } 
    
    /* A recursive function to insert a new key in BST */
    TrieNode insert(int[] populationId, double vol) { 
        int n = Global.population.length;
        int max_key = n - 1;
        int index;
        int level;
 
        /* If the tree is empty, return a new node */
        if (this.root == null) { 
            int key = -1;
            int keyes_required = populationId.length;
            int noChildren = max_key - key - keyes_required + 1;
            this.root = new TrieNode(key, noChildren); // -1 to denote it is root
        } 

        int min_possible = 0;
        
        TrieNode trieNode = this.root; 
        
        for (level = 0; level < populationId.length - 1; level++)  {
            int key = populationId[level];
            
            index = key - min_possible;
            if (trieNode.children[index] == null) { 
                int keyes_required = populationId.length - level - 1;
                int noChildren = max_key - key - keyes_required + 1;
                trieNode.children[index] = new TrieNode(key, noChildren);
            }
            trieNode = trieNode.children[index]; 
            min_possible = key+1;
        }
        
        int key = populationId[level];
        index = key - min_possible;
        trieNode.children[index] = new TrieNode(key, 1, vol, true);
        
        /*int key = populationId[level];
        trieNode.children[0] = new TrieNode(key, 1, vol, true);*/
        
        /* return the (unchanged) node pointer */
        return root; 
    } 
    
    
    public double search(int[] populationId) {
        int level; 
        int index; 
        TrieNode trieNode = this.root; 
        int min_possible = 0;

        if (trieNode == null) {
            return -1;
        }
        
        
        for (level = 0; level < populationId.length; level++) {
            int key = populationId[level];
            
            index = key - min_possible;
            
            if (trieNode.children[index] == null) {
                return -1;
            }
            
            min_possible = key+1;
            trieNode = trieNode.children[index]; 
        } 
       
        if (trieNode == null) {
            return -1;
        } else {
            return trieNode.volume; 
        }
    }
}

class TrieNode {
    int solutionId; 
    TrieNode[] children; 
    double volume;
    boolean isLeaf;
    
    public TrieNode() {
        this.isLeaf = false;
    }
    
    public TrieNode(int solutionId) { 
        this.solutionId = solutionId; 
        this.isLeaf = false;
    }   
    
    public TrieNode(int solutionId, int size) { 
        this.solutionId = solutionId; 
        this.children = new TrieNode[size];
        this.isLeaf = false;
    } 
    
    public TrieNode(int solutionId, int size, double volume, boolean isleaf) { 
        this.solutionId = solutionId; 
        this.children = new TrieNode[size];
        this.volume = volume; 
        this.isLeaf = true;
    } 

    public int getSolutionId() {
        return this.solutionId;
    }

    public TrieNode[] getChildren() {
        return this.children;
    }

    public double getVolume() {
        return volume;
    }

    public boolean isIsLeaf() {
        return this.isLeaf;
    }

    public void setSolutionId(int solutionId) {
        this.solutionId = solutionId;
    }

    public void setChildren(TrieNode[] children) {
        this.children = children;
    }

    public void setVolume(double volume) {
        this.volume = volume;
    }

    public void setIsLeaf(boolean isLeaf) {
        this.isLeaf = isLeaf;
    }

    @Override
    public String toString() {
        return "TrieNode{" + "solutionId=" + solutionId + ", children=" + children + ", volume=" + volume + ", isLeaf=" + isLeaf + '}';
    }
}
