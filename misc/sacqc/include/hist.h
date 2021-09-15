
typedef struct {
	float start;
	float stop;
	float center;
	double counts;
	double percent;
} Bin;

typedef struct {
	int nbin;
	float mode;
	float median;
	float mean;
	float xmax, xmax_percent;

	float maxstat;	/*** DEFAULT 5 the log10 of the mode, mean or median cannot be above this value ***/
	float minstat; 	/*** DEFAULT 0.5 the log10 value of the mode, mean, or median cannot be below this value ***/
	float max_x;	/*** DEFAULT 9 the max log10 amplitude cannot be above this value ***/
	float max_xper;	/*** DEFAULT 2 the max amp bin cannot contain more than 2 percent of the total amplitudes ***/
	char defmask;      /*** defining mask for the data 'Y' or 'N' ***/
	char reason[128];

	Bin *b;
} Histogram;

typedef struct {
        float ave, adev, svar, sdev;
        float skew, curt;
        float med, min, max;
}Statistics;

