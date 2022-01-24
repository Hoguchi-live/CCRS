#ifndef _CONFIG_H_
#define _CONFIG_H_

#define BASE_p "12037340738208845034383383978222801137092029451270197923071397735408251586669938291587857560356890516069961904754171956588530344066457839297755929645858769"
#define BASE_q "12037340738208845034383383978222801137092029451270197923071397735408251586669938291587857560356890516069961904754171956588530344066457839297755929645858769"

#define BASE_A	"10861338504649280383859950140772947007703646408372831934324660566888732797778932142488253565145603672591944602210571423767689240032829444439469242521864171"
#define BASE_B	"1"
#define BASE_t	"-147189550172528104900422131912266898599387555512924231762107728432541952979290"

#define NB_PRIMES 3
#define l_PRIMES_int int l_primes[NB_PRIMES] = {3, 5, 7};

/**
  Type for l-primes.
*/
typedef struct lprime_t{

	fmpz_t l;

	uint type; 		// Radical (0) or Velu (1)
	uint lbound, hbound;	// Bounds for the walk
	uint lorder, horder;	// Bounds on orders of the eigenvalues
	uint r;			//
	uint bkw;		// 1 if backward walking possible
} lprime_t ;

/// Config
typedef struct cfg_t{

	//// Base curve parameters
	const fq_ctx_t *F;		// base field
	const fmpz_t p;			// base field
	const fq_t A, B, t;		// base curve parameters

	//// l-primes parameters
	const uint nb_primes; 		// number of l-primes used
	const lprime_t[nb_primes];	// the l-primes
} cfg_t;

#endif

