/*
 * 2016 hankoor
 */

/**
Challenge description:
 *  Our marketing department has just negotiated a deal with several local merchants that will allow us to offer exclusive discounts on various products to our top customers every day. The catch is that we can only offer each product to one customer and we may only offer one product to each customer.

Each day we will get the list of products that are eligible for these special discounts.
We then have to decide which products to offer to which of our customers.
Fortunately, our team of highly skilled statisticians has developed an amazing mathematical model
for determining how likely a given customer is to buy an offered product
by calculating what we call the "suitability score" (SS).

The top-secret algorithm to calculate the SS between a customer and a product is this:

1. If the number of letters in the product's name is even 
	then the SS is the number of vowels (a, e, i, o, u, y) in the customer's name multiplied by 1.5.
2. If the number of letters in the product's name is odd 
	then the SS is the number of consonants in the customer's name.
3. If the number of letters in the product's name shares any common factors (besides 1) 
	with the number of letters in the customer's name then the SS is multiplied by 1.5.

Your task is to implement a program that assigns each customer a product to be offered
in a way that maximizes the combined total SS
across all of the chosen offers.
Note that there may be a different number of products and customers.
You may include code from external libraries as long as you cite the source.

Input sample:

Your program should accept as its only argument a path to a file.
Each line in this file is one test case.
Each test case will be a comma delimited set of customer names
	followed by a semicolon and then 
	a comma delimited set of product names. 
Assume the input file is ASCII encoded.
For example (NOTE: The example below has 3 test cases):

Jack Abraham,John Evans,Ted Dziuba;iPad 2 - 4-pack,Girl Scouts Thin Mints,Nerf Crossbow

Jeffery Lebowski,Walter Sobchak,Theodore Donald Kerabatsos,Peter Gibbons,Michael Bolton,Samir Nagheenanajar;Half & Half,Colt M1911A1,16lb bowling ball,Red Swingline Stapler,Printer paper,Vibe Magazine Subscriptions - 40 pack

Jareau Wade,Rob Eroh,Mahmoud Abdelkader,Wenyi Cai,Justin Van Winkle,Gabriel Sinkin,Aaron Adelson;Batman No. 1,Football - Official Size,Bass Amplifying Headphones,Elephant food - 1024 lbs,Three Wolf One Moon T-shirt,Dom Perignon 2000 Vintage

Output sample:

For each line of input, print out the maximum total score to two decimal places. For the example input above, the output should look like this:

21.00
83.50
71.25


https://www.codeeval.com/open_challenges/48/


 */
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class Main {
	final static String CONSONANTS = "bcdfghjklmnpqrstvwxz" + "BCDFGHJKLMNPQRSTVWXZ";
	final static String VOWELS = "AEIOUY" + "aeiouy";

	private static void loadFile(String fileName) {
		File file = new File(fileName);
		if (!file.canRead() || !file.isFile()) {
			System.exit(1);
		}
		BufferedReader in;
		String row = "";
		try {
			in = new BufferedReader(new FileReader(fileName));
			while ((row = in.readLine()) != null) {
				String[] names;
				String[] articles;
				int[] matching;
				if (row.split(";").length > 1) {
					names = row.split(";")[0].split(",");
					articles = row.split(";")[1].split(",");

					/*
					 * Please refer to the notes prior to nested class! It's a
					 * foreign class from the Massachusetts Institute of
					 * Technology I'm assuming, that the problem can be viewed
					 * as
					 * "maximum-weight matching in a complete bipartite graph"
					 */
					Main.MWBMatchingAlgorithm maxMatch = new Main.MWBMatchingAlgorithm(names.length, articles.length);
					/*
					 * 
					 */
					for (int i = 0; i < names.length; i++) {
						for (int j = 0; j < articles.length; j++) {
							double ss = 0;
							String name, article;
							name = names[i];
							article = articles[j];
							if (!numOfLettersIsOdd(article)) {
								ss = numOfVowels(name) * 1.5;
							} else {
								ss = numOfConsonants(name);
							}
							if (hasCommonFactor(numOfLetters(article), numOfLetters(name))) {
								ss *= 1.5;
							}
							maxMatch.setWeight(i, j, ss);
							}
					}
					matching = maxMatch.getMatching();
					System.out.printf("%.2f", maxMatch.calcWeightSumOfMatching());
					System.out.println();
				} else {
					System.out.println(00.00);
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	public static void main(String[] pathToFile) {
		if ((pathToFile != null) && (pathToFile.length == 1) && (!pathToFile[0].isEmpty())) {
			loadFile(pathToFile[0]);
		}
	}
	private static int numOfLetters(String str) {
		int numberOfLetters = 0;
		for (int i = 0; i < str.length(); i++) {
			if (Character.isLetter(str.charAt(i))) {
				numberOfLetters++;
			}
		}
		return numberOfLetters;
	}
	private static boolean numOfLettersIsOdd(String str) {
		if (numOfLetters(str) % 2 == 0) {
			return false;
		}
		return true;
	}
	private static int numOfVowels(String str) {
		int numberOfVowel = 0;
		for (int i = 0; i < str.length(); i++) {
			if (VOWELS.contains(str.charAt(i) + "")) {
				numberOfVowel++;
			}
		}
		return numberOfVowel;
	}
	private static int numOfConsonants(String str) {
		int numberOfConsonants = 0;
		for (int i = 0; i < str.length(); i++) {
			if (CONSONANTS.contains(str.charAt(i) + "")) {
				numberOfConsonants++;
			}
		}
		return numberOfConsonants;
	}
	private static boolean hasCommonFactor(int i, int y) {
		if ((i % 2) == 0 && (y % 2) == 0) {
			return true;
		}
		int max = Integer.max(i, y);
		int min = Integer.min(i, y);
		if ((min > 0) && (max % min == 0)) {
			return true;
		}
		if (min >= 3) {
			for (int j = 3; j <= min; j++) {
				if ((min % j == 0) && (max % j == 0)) {
					return true;
				}
			}
		}
		return false;
	}

	/*
	 * Copyright (c) 2007, Massachusetts Institute of Technology Copyright (c)
	 * 2005-2006, Regents of the University of California All rights reserved.
	 * 
	 * Redistribution and use in source and binary forms, with or without
	 * modification, are permitted provided that the following conditions are
	 * met:
	 *
	 * * Redistributions of source code must retain the above copyright notice,
	 * this list of conditions and the following disclaimer.
	 *
	 * * Redistributions in binary form must reproduce the above copyright
	 * notice, this list of conditions and the following disclaimer in the
	 * documentation and/or other materials provided with the distribution.
	 *
	 * * Neither the name of the University of California, Berkeley nor the
	 * names of its contributors may be used to endorse or promote products
	 * derived from this software without specific prior written permission.
	 *
	 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
	 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
	 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
	 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
	 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
	 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
	 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
	 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
	 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
	 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
	 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
	 */

	/**
	 * An engine for finding the maximum-weight matching in a complete bipartite
	 * graph. Suppose we have two sets <i>S</i> and <i>T</i>, both of size
	 * <i>n</i>. For each <i>i</i> in <i>S</i> and <i>j</i> in <i>T</i>, we have
	 * a weight <i>w <sub>ij</sub></i>. A perfect matching <i>X</i> is a subset
	 * of <i>S</i> x <i>T</i> such that each <i>i</i> in <i>S</i> occurs in
	 * exactly one element of <i>X</i>, and each <i>j</i> in <i>T</i> occurs in
	 * exactly one element of <i>X</i>. Thus, <i>X</i> can be thought of as a
	 * one-to-one function from <i>S</i> to <i>T</i>. The weight of <i>X</i> is
	 * the sum, over (<i>i</i>, <i>j</i>) in <i>X</i>, of <i>w<sub>ij</sub></i>.
	 * A BipartiteMatcher takes the number <i>n</i> and the weights <i>w
	 * <sub>ij</sub></i>, and finds a perfect matching of maximum weight.
	 *
	 * It uses the Hungarian algorithm of Kuhn (1955), as improved and presented
	 * by E. L. Lawler in his book <cite>Combinatorial Optimization: Networks
	 * and Matroids</cite> (Holt, Rinehart and Winston, 1976, p. 205-206). The
	 * running time is O(<i>n</i><sup>3</sup>). The weights can be any finite
	 * real numbers; Lawler's algorithm assumes positive weights, so if
	 * necessary we add a constant <i>c</i> to all the weights before running
	 * the algorithm. This increases the weight of every perfect matching by
	 * <i>nc</i>, which doesn't change which perfect matchings have maximum
	 * weight.
	 *
	 * If a weight is set to Double.NEGATIVE_INFINITY, then the algorithm will
	 * behave as if that edge were not in the graph. If all the edges incident
	 * on a given node have weight Double.NEGATIVE_INFINITY, then the final
	 * result will not be a perfect matching, and an exception will be thrown.
	 */
	static class MWBMatchingAlgorithm {
		/**
		 * Creates a BipartiteMatcher without specifying the graph size. Calling
		 * any other method before calling reset will yield an
		 * IllegalStateException.
		 */

		/**
		 * Tolerance for comparisons to zero, to account for floating-point
		 * imprecision. We consider a positive number to be essentially zero if
		 * it is strictly less than TOL.
		 */
		private static final double TOL = 1e-10;
		// Number of left side nodes
		int n;
		// Number of right side nodes
		int m;
		double[][] weights;
		double minWeight;
		double maxWeight;
		// If (i, j) is in the mapping, then sMatches[i] = j and tMatches[j] =
		// i.
		// If i is unmatched, then sMatches[i] = -1 (and likewise for tMatches).
		int[] sMatches;
		int[] tMatches;
		static final int NO_LABEL = -1;
		static final int EMPTY_LABEL = -2;
		int[] sLabels;
		int[] tLabels;
		double[] u;
		double[] v;
		double[] pi;
		List<Integer> eligibleS = new ArrayList<Integer>();
		List<Integer> eligibleT = new ArrayList<Integer>();

		public MWBMatchingAlgorithm() {
			n = -1;
			m = -1;
		}

		/**
		 * Creates a BipartiteMatcher and prepares it to run on an n x m graph.
		 * All the weights are initially set to 1.
		 */
		public MWBMatchingAlgorithm(int n, int m) {
			reset(n, m);
		}

		/**
		 * Resets the BipartiteMatcher to run on an n x m graph. The weights are
		 * all reset to 1.
		 */
		private void reset(int n, int m) {
			if (n < 0 || m < 0) {
				throw new IllegalArgumentException("Negative num nodes: " + n + " or " + m);
			}
			this.n = n;
			this.m = m;

			weights = new double[n][m];
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < m; j++) {
					weights[i][j] = 1;
				}
			}
			minWeight = 1;
			maxWeight = Double.NEGATIVE_INFINITY;

			sMatches = new int[n];
			tMatches = new int[m];
			sLabels = new int[n];
			tLabels = new int[m];
			u = new double[n];
			v = new double[m];
			pi = new double[m];

		}

		/**
		 * Sets the weight w<sub>ij</sub> to the given value w.
		 *
		 * @throws IllegalArgumentException
		 *             if i or j is outside the range [0, n).
		 */
		public void setWeight(int i, int j, double w) {
			if (n == -1 || m == -1) {
				throw new IllegalStateException("Graph size not specified.");
			}
			if ((i < 0) || (i >= n)) {
				throw new IllegalArgumentException("i-value out of range: " + i);
			}
			if ((j < 0) || (j >= m)) {
				throw new IllegalArgumentException("j-value out of range: " + j);
			}
			if (Double.isNaN(w)) {
				throw new IllegalArgumentException("Illegal weight: " + w);
			}

			weights[i][j] = w;
			if ((w > Double.NEGATIVE_INFINITY) && (w < minWeight)) {
				minWeight = w;
			}
			if (w > maxWeight) {
				maxWeight = w;
			}
		}

		/**
		 * Returns a maximum-weight perfect matching relative to the weights
		 * specified with setWeight. The matching is represented as an array arr
		 * of length n, where arr[i] = j if (i,j) is in the matching.
		 */
		public int[] getMatching() {
			if (n == -1 || m == -1) {
				throw new IllegalStateException("Graph size not specified.");
			}
			if (n == 0) {
				return new int[0];
			}
			ensurePositiveWeights();

			// Step 0: Initialization
			eligibleS.clear();
			eligibleT.clear();
			for (Integer i = 0; i < n; i++) {
				sMatches[i] = -1;

				u[i] = maxWeight; // ambiguous on p. 205 of Lawler, but see p.
									// 202

				// this is really first run of Step 1.0
				sLabels[i] = EMPTY_LABEL;
				eligibleS.add(i);
			}

			for (int j = 0; j < m; j++) {
				tMatches[j] = -1;

				v[j] = 0;
				pi[j] = Double.POSITIVE_INFINITY;

				// this is really first run of Step 1.0
				tLabels[j] = NO_LABEL;
			}

			while (true) {
				// Augment the matching until we can't augment any more given
				// the
				// current settings of the dual variables.
				while (true) {
					// Steps 1.1-1.4: Find an augmenting path
					int lastNode = findAugmentingPath();
					if (lastNode == -1) {
						break; // no augmenting path
					}

					// Step 2: Augmentation
					flipPath(lastNode);
					for (int i = 0; i < n; i++)
						sLabels[i] = NO_LABEL;

					for (int j = 0; j < m; j++) {
						pi[j] = Double.POSITIVE_INFINITY;
						tLabels[j] = NO_LABEL;
					}

					// This is Step 1.0
					eligibleS.clear();
					for (int i = 0; i < n; i++) {
						if (sMatches[i] == -1) {
							sLabels[i] = EMPTY_LABEL;
							eligibleS.add(new Integer(i));
						}
					}

					eligibleT.clear();
				}

				// Step 3: Change the dual variables

				// delta1 = min_i u[i]
				double delta1 = Double.POSITIVE_INFINITY;
				for (int i = 0; i < n; i++) {
					if (u[i] < delta1) {
						delta1 = u[i];
					}
				}

				// delta2 = min_{j : pi[j] > 0} pi[j]
				double delta2 = Double.POSITIVE_INFINITY;
				for (int j = 0; j < m; j++) {
					if ((pi[j] >= TOL) && (pi[j] < delta2)) {
						delta2 = pi[j];
					}
				}

				if (delta1 < delta2) {
					// In order to make another pi[j] equal 0, we'd need to
					// make some u[i] negative.
					break; // we have a maximum-weight matching
				}

				changeDualVars(delta2);
			}

			int[] matching = new int[n];
			for (int i = 0; i < n; i++) {
				matching[i] = sMatches[i];
			}
			return matching;
		}

		/**
		 * Tries to find an augmenting path containing only edges (i,j) for
		 * which u[i] + v[j] = weights[i][j]. If it succeeds, returns the index
		 * of the last node in the path. Otherwise, returns -1. In any case,
		 * updates the labels and pi values.
		 */
		int findAugmentingPath() {
			while ((!eligibleS.isEmpty()) || (!eligibleT.isEmpty())) {
				if (!eligibleS.isEmpty()) {
					int i = ((Integer) eligibleS.get(eligibleS.size() - 1)).intValue();
					eligibleS.remove(eligibleS.size() - 1);
					for (int j = 0; j < m; j++) {
						// If pi[j] has already been decreased essentially
						// to zero, then j is already labeled, and we
						// can't decrease pi[j] any more. Omitting the
						// pi[j] >= TOL check could lead us to relabel j
						// unnecessarily, since the diff we compute on the
						// next line may end up being less than pi[j] due
						// to floating point imprecision.
						if ((tMatches[j] != i) && (pi[j] >= TOL)) {
							double diff = u[i] + v[j] - weights[i][j];
							if (diff < pi[j]) {
								tLabels[j] = i;
								pi[j] = diff;
								if (pi[j] < TOL) {
									eligibleT.add(new Integer(j));
								}
							}
						}
					}
				} else {
					int j = ((Integer) eligibleT.get(eligibleT.size() - 1)).intValue();
					eligibleT.remove(eligibleT.size() - 1);
					if (tMatches[j] == -1) {
						return j; // we've found an augmenting path
					}

					int i = tMatches[j];
					sLabels[i] = j;
					eligibleS.add(new Integer(i)); // ok to add twice
				}
			}

			return -1;
		}

		/**
		 * Given an augmenting path ending at lastNode, "flips" the path. This
		 * means that an edge on the path is in the matching after the flip if
		 * and only if it was not in the matching before the flip. An augmenting
		 * path connects two unmatched nodes, so the result is still a matching.
		 */
		void flipPath(int lastNode) {
			while (lastNode != EMPTY_LABEL) {
				int parent = tLabels[lastNode];

				// Add (parent, lastNode) to matching. We don't need to
				// explicitly remove any edges from the matching because:
				// * We know at this point that there is no i such that
				// sMatches[i] = lastNode.
				// * Although there might be some j such that tMatches[j] =
				// parent, that j must be sLabels[parent], and will change
				// tMatches[j] in the next time through this loop.
				sMatches[parent] = lastNode;
				tMatches[lastNode] = parent;

				lastNode = sLabels[parent];
			}
		}

		void changeDualVars(double delta) {
			for (int i = 0; i < n; i++) {
				if (sLabels[i] != NO_LABEL) {
					u[i] -= delta;
				}
			}

			for (int j = 0; j < m; j++) {
				if (pi[j] < TOL) {
					v[j] += delta;
				} else if (tLabels[j] != NO_LABEL) {
					pi[j] -= delta;
					if (pi[j] < TOL) {
						eligibleT.add(new Integer(j));
					}
				}
			}
		}

		/**
		 * Ensures that all weights are either Double.NEGATIVE_INFINITY, or
		 * strictly greater than zero.
		 */
		private void ensurePositiveWeights() {
			// minWeight is the minimum non-infinite weight
			if (minWeight < TOL) {
				for (int i = 0; i < n; i++) {
					for (int j = 0; j < m; j++) {
						weights[i][j] = weights[i][j] - minWeight + 1;
					}
				}

				maxWeight = maxWeight - minWeight + 1;
				minWeight = 1;
			}
		}

		@SuppressWarnings("unused")
		private void printWeights() {
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < m; j++) {
					System.out.print(weights[i][j] + " ");
				}
				System.out.println("");
			}
		}

		public double calcWeightSumOfMatching() {
			int[] max = this.getMatching();
			double sum = 0;
			for (int i = 0; i < max.length; i++) {
				if (max[i] >= 0) {
					sum += weights[i][max[i]];
				}
			}
			return sum;
		}

	}

}
