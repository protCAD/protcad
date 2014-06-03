#define GLOBAL_DEBUG

#ifndef GLOBAL_DEBUG
	#define ASSERT(x)
#else
	#define ASSERT(x) \
	if (! (x)) \
	{ \
		cout << endl << "ERROR!! Assert " << #x << " failed" << endl; \
		cout << " on line " << __LINE__ << endl; \
		cout << " in file " << __FILE__ << endl; \
		char buffer(20); \
		cout << "continue? "; \
		cin >> buffer; \
	}
#endif
