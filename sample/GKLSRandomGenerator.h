#ifndef GKLSRandomGenerator_H
#define GKLSRandomGenerator_H

namespace gklsfunction
{
	namespace randomgenerator{
#define KK 100                     /* the long lag */
#define LL  37                     /* the short lag */
#define mod_sum(x,y) (((x)+(y))-(int)((x)+(y)))   /* (x+y) mod 1.0 */

#define TT  70   /* guaranteed separation between streams */
#define is_odd(s) ((s)&1)

#define QUALITY 1009 /* recommended quality level for high-res use */

#define NUM_RND 1009 /* size of the array of random numbers */

		class GKLSRandomGenerator final	{
		private:
			double* rnd_num; /* array of random numbers */
			double* ran_u;          /* the generator state */
		public:
			GKLSRandomGenerator();
			~GKLSRandomGenerator();

			void Initialize(long seed, double* rnd_num_mem, double* rand_condition);
			void GenerateNextNumbers();

			double GetRandomNumber(int indx) const; // indx must be less than NUM_RND
			double GetGeneratorState(int indx) const; // indx must be less than KK
		};
	}
}
#endif