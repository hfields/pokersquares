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
	private Random random = new Random(); // pseudorandom number generator for breaking ties when choosing moves
	private int[] plays = new int[NUM_POS]; // positions of plays so far (index 0 through numPlays - 1) recorded as integers using row-major indices.
	// row-major indices: play (r, c) is recorded as a single integer r * SIZE + c (See http://en.wikipedia.org/wiki/Row-major_order)
	// From plays index [numPlays] onward, we maintain a list of yet unplayed positions.
	private int numPlays = 0; // number of Cards played into the grid so far
	private int recursions = 0; // debugging variable for counting number of recursive calls of expectimax
	private PokerSquaresPointSystem system; // point system
	private int depthLimit = 2; // default depth limit for expectimax search
	private Card[][] grid = new Card[SIZE][SIZE]; // grid with Card objects or null (for empty positions)
	private Card[] simDeck = Card.getAllCards(); // a list of all Cards. As we learn the index of cards in the play deck,
	                                             // we swap each dealt card to its correct index.  Thus, from index numPlays 
												 // onward, we maintain a list of undealt cards for expectimax calculations.
	private int[] legalPlayList = new int[NUM_POS]; // stores legal play lists indexed by numPlays (depth)
	
	/**
	 * Create an expectimax player that performs depth-limited heuristic expectimax
	 * search to a depth of 3.
	 */
	public ExpectimaxPlayer() {
	}
	
	/**
	 * Create an expectimax player that performs depth-limited heuristic expectimax
	 * search to a given depth limit.
	 * @param depthLimit depth limit for expectimax search
	 */
	public ExpectimaxPlayer(int depthLimit) {
		this.depthLimit = depthLimit;
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
			
			for (int i = 0; i < remainingPlays; i++) { // for each legal play position
				int play = legalPlayList[i];		
				double expectimax = heuristicExpectimax(card, play / SIZE, play % SIZE, this.depthLimit);  // find expectimax value of current play
				
				// Update max expectimax value/best plays list, if necessary
				if (expectimax >= maxExpectimax) {
					if (expectimax > maxExpectimax) {
						bestPlays.clear();
					}
					
					bestPlays.add(play);
					maxExpectimax = expectimax;
				}
		
			}
			int bestPlay = bestPlays.get(random.nextInt(bestPlays.size())); // choose a best play (breaking ties randomly)
			// update our list of plays, recording the chosen play in its sequential position; all onward from numPlays are empty positions
			int bestPlayIndex = numPlays;
			while (plays[bestPlayIndex] != bestPlay)
				bestPlayIndex++;
			plays[bestPlayIndex] = plays[numPlays];
			plays[numPlays] = bestPlay;
		}
		
		// Print number of calls to expectimax required for the move
		System.out.print("Calls to expectimax: ");
		System.out.println(recursions);
		
		int[] playPos = {plays[numPlays] / SIZE, plays[numPlays] % SIZE}; // decode it into row and column
		makePlay(card, playPos[0], playPos[1]); // make the chosen play (not undoing this time)
		return playPos; // return the chosen play
	}
	
	/**
	 * Returns a heuristic evaluation of the given game state
	 */
	private double heuristic() {
		double eval = 0;
		
		int[] handScores = system.getHandScores(grid);
		
		for (int i = 0; i < 2 * SIZE; i++) {
			eval += handScores[i];
		}
		
		return eval;
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
			System.arraycopy(plays, numPlays, legalPlayList, 0, remainingPlays);
			
			for (int i = 0; i < remainingPlays; i++) {
				double newExpectimax = 0;
				int newPlay = legalPlayList[i];
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
	 * Demonstrate ExpectimaxPlayer play with Ameritish point system.
	 * @param args (not used)
	 */
	public static void main(String[] args) {
		PokerSquaresPointSystem system = PokerSquaresPointSystem.getBritishPointSystem();
		System.out.println(system);
		new PokerSquares(new ExpectimaxPlayer(), system).play(); // play a single game
	}

}