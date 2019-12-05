import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;
import java.util.Stack;

/**
 * ExpectimaxPlayer - an implementation of the player interface for PokerSquares that 
 * makes placements based on depth-limited heuristic expectimax 
 * Authors: Huey Fields and Sol Cruz, with inspiration from Erin Talvitie
 * Code modified from GreedyMCPlayer implementation by Todd W. Neller
 */
public class ExpectimaxPlayer implements PokerSquaresPlayer {
	private final int SIZE = 5; // number of rows/columns in square grid
	private final int NUM_POS = SIZE * SIZE; // number of positions in square grid
	private final int NUM_CARDS = Card.NUM_CARDS; // number of cards in deck
	private boolean debug = true; // boolean determining whether or not debug statements will print
	private Random random = new Random(); // pseudorandom number generator for breaking ties when choosing moves
	private int[] plays = new int[NUM_POS]; // positions of plays so far (index 0 through numPlays - 1) recorded as integers using row-major indices.
	// row-major indices: play (r, c) is recorded as a single integer r * SIZE + c (See http://en.wikipedia.org/wiki/Row-major_order)
	// From plays index [numPlays] onward, we maintain a list of yet unplayed positions.
	private int numPlays = 0; // number of Cards played into the grid so far
	private int recursions = 0; // debugging variable for counting number of recursive calls of expectimax
	private PokerSquaresPointSystem system; // point system
	private double recursionLimit = 1000000; // rough limit for number of calls that should be made to expectimax to find a single move.
										  // enforces a changing maximum depth as the game progresses
	private Card[][] grid = new Card[SIZE][SIZE]; // grid with Card objects or null (for empty positions)
	private Card[] simDeck = Card.getAllCards(); // a list of all Cards. As we learn the index of cards in the play deck,
	                                             // we swap each dealt card to its correct index.  Thus, from index numPlays 
												 // onward, we maintain a list of undealt cards for expectimax calculations.
	private int[] legalPlayList = new int[NUM_POS]; // stores legal play lists indexed by numPlays (depth)
	
	/**
	 * Create an expectimax player that performs depth-limited heuristic expectimax
	 * search with a rough recursion limit of 1000000 expectimax calls.
	 */
	public ExpectimaxPlayer() {
	}
	
	/**
	 * Create an expectimax player that performs depth-limited heuristic expectimax
	 * search with depth determined by the given recursion limit.
	 * @param recursionLimit recursion limit to inform depth limit of expectimax search
	 */
	public ExpectimaxPlayer(int recursionLimit, boolean debug) {
		this.recursionLimit = recursionLimit;
		this.debug = debug;
	}
	
	/* (non-Javadoc)
	 * @see PokerSquaresPlayer#setPointSystem(PokerSquaresPointSystem, long)
	 */
	@Override
	public void setPointSystem(PokerSquaresPointSystem system, long millis) {
		// The expectimax player is meant to work exclusively with the UK scoring system
		this.system = system;
	}
	
	/* (non-Javadoc)
	 * @see PokerSquaresPlayer#init()
	 */
	@Override
	public void init() {
		// clear grid
		for (int row = 0; row < SIZE; row++)
			for (int col = 0; col < SIZE; col++)
				grid[row][col] = null;
		// reset numPlays
		numPlays = 0;
		// (re)initialize list of play positions (row-major ordering)
		for (int i = 0; i < NUM_POS; i++)
			plays[i] = i;
	}

	/* (non-Javadoc)
	 * @see PokerSquaresPlayer#getPlay(Card, long)
	 */
	@Override
	public int[] getPlay(Card card, long millisRemaining) {
		/*
		 * With this algorithm, the player chooses the legal play that has the highest expected score outcome.
		 * This outcome is estimated as follows:
		 *   For each move, we calculate the expectimax value of the move to the given depth limit.
		 *   If a terminal state is reached before    
		 */
		
		// match simDeck to actual play event; in this way, all indices forward from the card contain a list of 
		//   undealt Cards in some permutation.
		int cardIndex = numPlays;
		while (!card.equals(simDeck[cardIndex]))
			cardIndex++;
		simDeck[cardIndex] = simDeck[numPlays];
		simDeck[numPlays] = card;
		
		if (numPlays < 24) { // not the forced last play
			// compute average time per move evaluation
			int remainingPlays = NUM_POS - numPlays; // ignores triviality of last play to keep a conservative margin for game completion
			// copy the play positions (row-major indices) that are empty
			System.arraycopy(plays, numPlays, legalPlayList, 0, remainingPlays);
			ArrayList<Integer> bestPlays = new ArrayList<Integer>(); // all plays yielding the maximum average score 
			double maxExpectimax = Double.NEGATIVE_INFINITY;
			
			// Compute a rough estimate of recursions needed to calculate expectimax values for this move
			// Use this to compute depth limit for expectimax
			double recursionEstimate = remainingPlays;
			int depthLimit = 0;
			for (int depth = 0; depth < remainingPlays; depth++) {
				
				if (debug) {
					System.out.println(recursionEstimate);
				}
				
				int remPlays = remainingPlays - (depth + 1);
				int remCards = NUM_CARDS - (NUM_POS - remPlays);
				recursionEstimate = recursionEstimate + recursionEstimate * remPlays * remCards;
				
				// Set depthLimit and break if we exceed recursion limit
				if (recursionEstimate > this.recursionLimit) {
					depthLimit = depth;
					break;
				}
				
				// Or set depthLimit to max depth if we reach the end
				else if (depth == remainingPlays - 1) {
					depthLimit = depth;
					break;
				}
			}
			
			if (debug) {
				System.out.print("Depth limit at move ");
				System.out.print(numPlays + 1);
				System.out.print(": ");
				System.out.println(depthLimit);
				System.out.print("Number of recursions for 1 deeper depth limit: ");
				System.out.println(recursionEstimate);
			}
			
			for (int i = 0; i < remainingPlays; i++) { // for each legal play position
				int play = legalPlayList[i];
				double expectimax = heuristicExpectimax(card, play / SIZE, play % SIZE, depthLimit);  // find expectimax value of current play
				
				// Update max expectimax value/best plays list, if necessary
				if (expectimax >= maxExpectimax) {
					if (expectimax > maxExpectimax) {
						bestPlays.clear();
					}
					
					bestPlays.add(play);
					maxExpectimax = expectimax;
				}
		
			}
			
			if (debug) {
				System.out.println(maxExpectimax);
				System.out.println(bestPlays);
			}
			
			int bestPlay = bestPlays.get(random.nextInt(bestPlays.size())); // choose a best play (breaking ties randomly)
			// update our list of plays, recording the chosen play in its sequential position; all onward from numPlays are empty positions
			int bestPlayIndex = numPlays;
			while (plays[bestPlayIndex] != bestPlay)
				bestPlayIndex++;
			plays[bestPlayIndex] = plays[numPlays];
			plays[numPlays] = bestPlay;
		}
		
		if (debug) {
			// Print number of calls to expectimax required for the move
			System.out.print("Calls to expectimax: ");
			System.out.println(recursions);
		}
		
		int[] playPos = {plays[numPlays] / SIZE, plays[numPlays] % SIZE}; // decode it into row and column
		makePlay(card, playPos[0], playPos[1]); // make the chosen play (not undoing this time)
		return playPos; // return the chosen play
	}
	
	
	/**
	 * Calculates the expectimax value of a given placement of a card on the grid. 
	 * @param card the card being placed on the board
	 * @param row the row that the card is being placed in
	 * @param col the column that the card is being placed in
	 * @param depthLimit recursion limit. If this limit is exceeded, a heuristic value
	 * 	will be calculated and returned as the expectimax value.
	 */
	private double heuristicExpectimax(Card card, int row, int col, int depthLimit) {
		recursions++;
		
		double expectimax = Double.NEGATIVE_INFINITY;
		
		// Temporarily make the play indicated by card, row, and col
		makePlay(card, row, col);
		
		
		// If we are at the max depth of recursion, find value with heuristic evaluation function
		if (depthLimit == 0) {
			expectimax = heuristic();
			
			// Undo temporary move
			undoPlay();
			return expectimax;
		}
		
		else {
			int remainingPlays = NUM_POS - numPlays;
			
			// If the state is terminal, simply return score of current grid
			if (remainingPlays == 0) {
				expectimax = system.getScore(grid);
				
				// Undo temporary move
				undoPlay();
				return expectimax;
			}
			
			// copy the play positions (row-major indices) that are empty
			int[] legalPlays = new int[NUM_POS];
			System.arraycopy(plays, numPlays, legalPlays, 0, remainingPlays);
			
			for (int i = 0; i < remainingPlays; i++) {
				double newExpectimax = 0;
				int newPlay = legalPlays[i];
				int newRow = newPlay / 5;
				int newCol = newPlay % 5;
				
				// Sum expectimax values for each possible next card
				for (int j = numPlays; j < NUM_CARDS; j++) {
					newExpectimax += heuristicExpectimax(simDeck[j], newRow, newCol, depthLimit - 1);
				}
				
				// Divide by number of cards remaining to find average value (since there is an equal chance of getting each card)
				newExpectimax = newExpectimax / (NUM_CARDS - numPlays);
				
				// Update final expectimax score if necessary
				if (newExpectimax > expectimax) {
					expectimax = newExpectimax;
				}
			}
			
			// Undo temporary move
			undoPlay();
			return expectimax;
		}
	}
	
	/**
	 * Returns a heuristic evaluation of the given game state
	 */
	private double heuristic() {
		double eval = 0;
		
		int[] handScores = system.getHandScores(grid);
		
		for (int i = 0; i < 2 * SIZE; i++) {
			int handScore = handScores[i];
			
			// If a hand is already made, it contributes its score to the heuristic val
			if (handScore > 0) {
				eval += handScores[i];
			}
			
			// Otherwise, evaluate potential "goodness" of yet to be made hands
			else {
				// Check rows
				if (i < 5) {	
					//eval += 0;
					eval += rowEval(i);
				}
				
				// Check columns
				else {
					eval += colEval(i-5);
				}
			}
		}
		
		return eval;
	}
	
	/**
	 * Calculates the score that a given row will contribute to the heuristic value. 
	 * @param row the row being evaluated
	 */
	private double rowEval(int row) {
		double score = 0;
		
		// Add cards in the row to ArrayList of cards
		ArrayList<Card> cards = new ArrayList<Card>();
		
		for (int col = 0; col < SIZE; col++) {
			Card card = grid[row][col];
			
			if (card != null) {
				cards.add(grid[row][col]);
			}
		}
		
		// Calculate straight draw potential
		//score = getFlushDrawScore(cards);
		
		return score;
	}
	
	/**
	 * Calculates the score that a given column will contribute to the heuristic value. 
	 * @param col the column being evaluated
	 */
	private double colEval(int col) {
		double score = 0;
		
		// Add cards in the column to ArrayList of cards
		ArrayList<Card> cards = new ArrayList<Card>();
		
		for (int row = 0; row < SIZE; row++) {
			Card card = grid[row][col];
			
			if (card != null) {
				cards.add(grid[row][col]);
			}
		}
		
		// Calculate straight draw potential
		score += getStraightDrawScore(cards);
		
		return score;
	}
	
	/**
	 * Evaluates the flush draw potential of a given ArrayList of Cards.
	 * Every card that contributes to the flush draw adds 1 to the final
	 * result. A missed flush draw (cards of two different suits) returns 0. 
	 * @param cards the cards to evaluate the flush draw potential of
	 */
	private double getFlushDrawScore(ArrayList<Card> cards) {
		int[] scores = system.getScoreTable();
		int flushScore = scores[PokerHand.FLUSH.id];
		
		// Initialize variables to check for flush draws
		int currentSuit = -1; //-1 indicates no suit set yet
		int suitCount = 0;
		double score = 0;
		
		for (Card card : cards) {
			// Update flush draw if this is the first card in row or it matches the last card's suit
			if (card.getSuit() == currentSuit || currentSuit == -1) {
				suitCount += 1;
				currentSuit = card.getSuit();
			}
			
			// If suit does not match, flush draw is no longer possible
			else {
				suitCount = 0;
				break;
			}
		}
		
		// If we have more than one card of the same suit, add appropriate
		// fraction of full flush score
		if (suitCount > 1) {
			score = flushScore * ((float)suitCount / SIZE);
		}
		
		// If we have one card of the same suit, set to half of appropriate
		// fraction of full flush score. This is to encourage agent not to
		// ruin flush draws and continually stack up same-suit cards on the
		// same row
		else if (suitCount == 1) {
			score = flushScore * ((float)suitCount / SIZE / 2);
		}
		
		return score;
	}
	
	/**
	 * Evaluates the straight draw potential of a given ArrayList of Cards.
	 * Every card that contributes to the straight draw adds to the final result. 
	 * A missed straight draw (cards of two different suits) returns 0. 
	 * @param cards the cards to evaluate the flush draw potential of
	 */
	private double getStraightDrawScore(ArrayList<Card> cards) {
		int[] scores = system.getScoreTable();
		int straightScore = scores[PokerHand.STRAIGHT.id];
		
		// Initialize variables to check for straight draws
		int cardCount = 0;
		int straightDraw = 0;
		int maxRank = -1;
		int minRank = 13;
		int lowRank = 0;
		int highRank = 12;
		ArrayList<Integer> usedRanks = new ArrayList<Integer>();
		double score = 0;
		
		for (Card card : cards) {
			int rank = card.getRank();
			
			// If there is only an ace in the straight draw, allow it to count as
			// high or low, depending on the new card
			if (maxRank == 0 && minRank == 0) {
				// If card is 10 or higher, count ace as high
				if (rank >= 9) {
					maxRank = 13;
					minRank = 13;
				}
				
				// If card is 5 or lower, count ace as low
				else if (rank <= 5) {
					maxRank = 0;
					minRank = 0;
				}
				
				// Otherwise, straight draw is missed
				else {
					straightDraw = 0;
					break;
				}
			}
			
			// If there is at least a 10 in the straight draw, treat ace as high
			if (maxRank >= 9 && rank == 0) {
				rank = 13;
			}
			
			// Update straight draw if this is the first card in col or it contributes to the straight
			if (straightDraw == 0 || (rank >= lowRank && rank <= highRank && !usedRanks.contains(rank))) {
				straightDraw += 1;
				
				// Update max and min ranks in straight draw, if needed
				if (rank > maxRank) {
					maxRank = rank;
				}
				
				if (rank < minRank) {
					minRank = rank;
				}
				
				// Calculate new bounds for straight
				lowRank = minRank - (SIZE - 1 - (maxRank - minRank));
				highRank = maxRank + (SIZE - 1 - (maxRank - minRank));
				
				// Ensure that bounds do not exceed limits
				if (lowRank < 0) {
					lowRank = 0;
				}
				
				if (highRank > 13) {
					highRank = 13;
				}
			}
			
			// If straight draw is missed, break immediately
			else {
				straightDraw = 0;
				break;
			}
		}
		
		// If we have more than one card in a straight, add appropriate
		// fraction of full straight score
		if (straightDraw > 1) {
			score = straightScore * ((float)straightDraw / SIZE);
		}
		
		// If we have one card of the same suit, set to half of appropriate
		// fraction of full straight score. This is to encourage agent not to
		// ruin straight draws and continually stack up same-suit cards on the
		// same row
		else if (straightDraw == 1) {
			score = straightScore * ((float)straightDraw / SIZE / 2);
		}
		
		return score;
	}
	
	public void makePlay(Card card, int row, int col) {
		// match simDeck to event
		int cardIndex = numPlays;
		while (!card.equals(simDeck[cardIndex]))
			cardIndex++;
		simDeck[cardIndex] = simDeck[numPlays];
		simDeck[numPlays] = card;
		
		// update plays to reflect chosen play in sequence
		grid[row][col] = card;
		int play = row * SIZE + col;
		int j = 0;
		while (plays[j] != play)
			j++;
		plays[j] = plays[numPlays];
		plays[numPlays] = play;
		
		// increment the number of plays taken
		numPlays++;
	}

	public void undoPlay() { // undo the previous play
		numPlays--;
		int play = plays[numPlays];
		grid[play / SIZE][play % SIZE] = null;
	}

	/* (non-Javadoc)
	 * @see PokerSquaresPlayer#getName()
	 */
	@Override
	public String getName() {
		return "ExpectimaxPlayer";
	}

	/**
	 * Demonstrate ExpectimaxPlayer play with British point system.
	 * @param args (not used)
	 */
	public static void main(String[] args) {
		PokerSquaresPointSystem system = PokerSquaresPointSystem.getBritishPointSystem();
		System.out.println(system);
		new PokerSquares(new ExpectimaxPlayer(), system).play(); // play a single game
	}

}