#ifndef LOG_TOOL_CLASSES
#define LOG_TOOL_CLASSES


struct Date
	{	
		int day; 
		int month; 
		int year; 
	
		Date& operator= (const Date& date)
		{
			day = date.day;
			month = date.month;
			year = date.year;
			return *this;
		}
	};


	class BaseData
	{
	public:
		int ver_num;
		int rel_num;
		int log_num;
		Date date;
		mystr time;
		mystr user;									
		bool public_internal;				//defines whether a log/bug is public or internal
		mystr entry;								//the log entry / description

			
		BaseData() : ver_num(0), rel_num(0), log_num(0), date(), time(), user(), public_internal(false), entry() 
		{}

		BaseData(BaseData& data)
		{
			CopyFrom(data);
		}

		void CopyFrom(const BaseData& data)
		{
			ver_num = data.ver_num;
			rel_num = data.rel_num;
			log_num = data.log_num;
			date = data.date;
			time = data.time;
			user = data.user;
			public_internal = data.public_internal;
			entry = data.entry;
		}
	};
	


	class LogData : public BaseData
	{
	public:
		mystr type;									//the log type ... new feature (f), change (c), extension (e)
		
	
		LogData() : BaseData(), type()
		{ }

		LogData(LogData& data)
		{
			CopyFrom(data);
		}

		void CopyFrom(const LogData& data)
		{
			BaseData::CopyFrom(data);
			type = data.type;
		}
	};

	////struct defining the data of the log
	//struct LogData 
	//{
	//	int ver_num;
	//	int rel_num;
	//	int log_num;
	//	Date date;
	//	mystr time;
	//	mystr user;									
	//	mystr type;									//the log type ... new feature (f), change (c), extension (e)
	//	mystr entry;								//the log entry / description
	//	bool public_internal;				//defines whether a log is public or internal
	//};

	class BugData: public BaseData
	{
	public:
		mystr solution;							//possibly the solution to the bug....
		bool solved;								//bug solved??
	
		BugData() : BaseData()
		{ }

		BugData(BugData& data)
		{
			CopyFrom(data);
		}

		void CopyFrom(const BugData& data)
		{
			BaseData::CopyFrom(data);
			solution = data.solution;
			solved = data.solved;
		}
	};

	//struct BugData 
	//{
	//	int ver_num;
	//	int rel_num;
	//	int bug_num;
	//	Date date;
	//	mystr time;
	//	mystr user;									
	//	mystr entry;								//the log entry / description
	//	mystr solution;							//possibly the solution to the bug....
	//	bool solved;								//bug solved??	
	//};


	//Sorts an array a(1..a.Length()) into ascending numerical order by Shell’s method (diminishing increment
//sort). a is replaced on output by its sorted rearrangement.

template <typename T, typename T2>
void Quicksort(TArray<T>& a, TArray<T2>& b)
{
	int n = a.Length();
	int i,j,inc;
	T v;
	T2 vb;
	inc=1; //Determine the starting increment.
	do 
	{
		inc *= 3;
		inc++;
	} while (inc <= n);
	do 
	{ //Loop over the partial sorts.
		inc /= 3;
		for (i=inc;i<n;i++) 
		{ //Outer loop of straight insertion.
			v = a.Get0(i);
			vb = b.Get0(i);
			j=i;
			while (a.Get0(j-inc) > v) 
			{ //Inner loop of straight insertion.
				a.Elem0(j) = a.Get0(j-inc);
				b.Elem0(j) = b.Get0(j-inc);
				j -= inc;
				if (j < inc) break;
			}
			a.Elem0(j)=v;
			b.Elem0(j)=vb;
		}
	} while (inc > 1);

}

#endif