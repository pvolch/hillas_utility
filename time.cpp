/*time.cpp */
using namespace std;

class time_cam
{
public:
int year;
int month;
int day;
int hour;
int minute;
int second;
int mlsec;
int mksec;
int nsec;

time_cam(int Ut, int y, int mon, int d, int h, int min, int sec, int mls, int mks, int nan){
	UTS = Ut;
	year = y;
	month = mon;
	day = d;
	hour = h;
	minute = min;
	second = sec;
	mlsec = mls;
	mksec = mks;
	nsec = nan;
}
time_cam(){
}
double get_unix_time() {
	time_t rawtime;
	double nana = 0;
	struct tm * timeinfo;
	/* get current timeinfo: */
	time ( &rawtime ); //or: rawtime = time(0);
	/* convert to struct: */
	timeinfo = localtime ( &rawtime );
	//cout << mktime (timeinfo) << endl;
	/* now modify the timeinfo to the given date: */
	timeinfo->tm_year   = year - 1900;//year = 20, 2020 - 1900
	timeinfo->tm_mon    = month - 1;//months since January - [0,11]
	timeinfo->tm_mday   = day;      //day of the month - [1,31]
	timeinfo->tm_hour   = hour+UTS;     //hours since midnight - [0,23]
	timeinfo->tm_min    = minute;      //minutes after the hour - [0,59]
	timeinfo->tm_sec    = second;      //seconds after the minute - [0,59]
	if(nsec >= 500) {
		nana = 1;
	}
	return mktime(timeinfo) + double(mlsec) * 1e-3 + double(mksec + nana) * 1e-6;
}

unsigned int get_nsec(){
	return mlsec*1e6 + mksec*1e3 + nsec;
}

void set_time(string date, string tim){
	int dat[6];
	string source;
	char delim[] = {':',':',',','.','.','.'};
	day = stoi(date.substr(0,2)); //day
	month = stoi(date.substr(2,2)); //months
	year = 2000 + stoi(date.substr(4,2)); //year

	stringstream iss(tim);
	for(int i = 0; i < 6; i++) {
		getline(iss, source, delim[i]);
		istringstream iss(source);
		iss >> dat[i];
		//cout << dat[i] << ":";
	}
	hour = dat[0];
	minute = dat[1];
	second = dat[2];
	mlsec = dat[3];
	mksec = dat[4];
	nsec = dat[5];
	//  cout << endl;
}

void get_human_string_data(){
	printf("%02d.%02d.%02d\n", day, month, year);
}

void get_human_string_time(){
	printf("%02d:%02d:%02d,%03d.%03d.%03d\n", hour, minute, second, mlsec, mksec, nsec);
}

void get_human_string_data_time(){
	printf("%02d.%02d.%02d\t%02d:%02d:%02d,%03d.%03d.%03d\n", day, month, year, hour, minute, second, mlsec, mksec, nsec);
}

private:
int UTS = 3;
};
