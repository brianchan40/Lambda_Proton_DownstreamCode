#include <my_func.cpp>

void optimizing_s_b_ratio()
{
    double bmin, bmax;
    double bmin_curr = 1.11;
    double bmax_curr = 1.117;
    double significance = 0;
    double significance_curr = 0;

    for(int i = 0; i < 10; i++)
    {
        for(int j = 0; j < 10; j++)
        {
        	cout << "bmin_curr = " << bmin_curr << ", bmax_curr = " << bmax_curr << endl;
            significance_curr = my_func(bmin_curr, bmax_curr);

            if(significance_curr > significance)
            {
                bmin = bmin_curr;
                bmax = bmax_curr;
                significance = significance_curr;
            }

            bmin_curr += 0.0005;
        }

        bmin_curr = 1.11;
        bmax_curr += 0.0005;
    }

    cout << "bmin = " << bmin << ", bmax = " << bmax << endl;
    cout << "significance = " << significance << endl;
}