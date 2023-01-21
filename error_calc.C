using namespace std;

double error_mult(double a, double b, double a_err, double b_err){
	double f_abs = fabs(a*b);
	double a_err_pow = pow((a_err/a), 2);
	double b_err_pow = pow((b_err/b), 2);

	return f_abs * sqrt(a_err_pow + b_err_pow);
}

double error_divide(double a, double b, double a_err, double b_err){
	double f_abs = fabs(a/b);
	double a_err_pow = pow((a_err/a), 2);
	double b_err_pow = pow((b_err/b), 2);

	return f_abs * sqrt(a_err_pow + b_err_pow);
}

double error_add(double a_err, double b_err){
	double a_err_pow = pow(a_err, 2);
	double b_err_pow = pow(b_err, 2);

	return sqrt(a_err_pow + b_err_pow);
}

double purity_correction(double signal, double background, double purity_stat)
{
    return ((signal - (1.0 - purity_stat) * background) / purity_stat);
}

double purity_correction_error(double signal_error, double background_error, double purity_stat)
{
    return (sqrt(pow((1.0 / purity_stat), 2) * pow(signal_error, 2) + pow(((1 - purity_stat) / purity_stat), 2) * pow(background_error, 2)));
}

double efficiency_sum_pt(double signal, double number, string lam, int cent_id, int pt_id)
{
    if(lam == "lam") return signal * 1.0 / (double)lam_eff[cent_id][pt_id + 5] * number;
    else if(lam == "antilam") return signal * 1.0 / (double)antilam_eff[cent_id][pt_id + 5] * number;
}

double efficiency_sum_pt_error(double old_error, double signal, double signal_error, double number, string lam, int cent_id, int pt_id)
{
    double temp112_ss_error2 = 0;
    if(lam == "lam") temp112_ss_error2 = fabs(signal * 1.0 / (double)lam_eff[cent_id][pt_id + 5] * number) * sqrt(pow((signal_error / signal), 2) + 1.0 / number);
    else if(lam == "antilam") temp112_ss_error2 = fabs(signal * 1.0 / (double)antilam_eff[cent_id][pt_id + 5] * number) * sqrt(pow((signal_error / signal), 2) + 1.0 / number);

    return sqrt(pow(old_error, 2) + pow(temp112_ss_error2, 2));
}

double simple_efficiency_sum_pt_error(double signal_error, string lam, int cent_id, int pt_id){
	double temp112_ss_error2 = 0;

    if(lam == "lam") temp112_ss_error2 = signal_error/(double)lam_eff[cent_id][pt_id + 5];
    else if(lam == "antilam") temp112_ss_error2 = signal_error/(double)antilam_eff[cent_id][pt_id + 5];

    return temp112_ss_error2;
}

double number_sum_pt(double number, string lam, int cent_id, int pt_id)
{
    if(lam == "lam") return 1.0 / (double)lam_eff[cent_id][pt_id + 5] * number;
    else if(lam == "antilam") return 1.0 / (double)antilam_eff[cent_id][pt_id + 5] * number;
}

double number_sum_pt_error(double old_error, double number, string lam, int cent_id, int pt_id)
{
    if(lam == "lam") return sqrt(pow(old_error, 2) + pow(1.0 / (double)lam_eff[cent_id][pt_id + 5] * sqrt(number), 2));
    else if(lam == "antilam") return sqrt(pow(old_error, 2) + pow(1.0 / (double)antilam_eff[cent_id][pt_id + 5] * sqrt(number), 2));
}