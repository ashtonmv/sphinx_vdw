#include <SxIterator.h>
#include <SxArray.h>
#include <SxList.h>
#include <SxRange.h>
#include <SxError.h>
#include <SxString.h>
#include <SxCLI.h>
#include <SxLog.h>
void foo2() { int i = 2; SX_CHECK(i == 1, i, i, i, i); }

void foo1() { foo2(); }

template<class Fn>
void myDbgMsg(Fn f) { f(int(0)); }

#define MYDBGMSG(msg) myDbgMsg([](auto) { cout <<"DBG: "<<__LINE__<<" " <<msg; });

int main(int argc, char **argv)
{
	SxCLI cli(argc, argv);
	cli.finalize();
	//wchar_t *a = new wchar_t[5];
	//SxString b("asdqesda");
	//a[0] = 'a'; a[1] = 'b'; a[2] = 'c'; a[3] = 'd'; a[4] = '\0';
	//const wchar_t *aa = L"";

	//char *a = u8"παράδειγμα";
	//wchar_t *b = L"  παράδειγμα";
	//SxString s("  hello string");
    //std::wcout << a << b << s << std::endl;

	//char *a = L"παράδειγμα";

	//SX_DBG_MSG(a);

	SX_DBG_MSG(u8"παράδειγμα");
	SX_LOG(u8"παράδειγμα");
	SX_LOG(u8"пример");
	SX_LOG(u8"例");
	SX_LOG(u8"ตัวอย่าง");

	//std::wcerr << L"παράδειγμα" << std::endl;
	//std::wcerr << L"пример" << std::endl;

	/*SX_DBG_MSG(L"παράδειγμα");
	SX_LOG(L"παράδειγμα");
	SX_LOG(L"пример");
	SX_LOG(L"例");
	SX_LOG(L"ตัวอย่าง");*/

	/*SxList<int> list;
	SX_DBG_MSG ("hi",123);
	list << 1 << 2 << 3;

	SX_BREAK;

	//   myDbgMsg ([](auto) { cout<<"DBG: "<<__FILE__<<":"<<__LINE__<<endl;});
	cout << "RES=" << *(list.begin().SX_BREAK
	.SX_DBG_MSG(_it, "aaa: " << *_it)
	.SX_BREAK
	.print()) << endl;*/

	//   SxList<int>::Iterator iter;
	//   if (list.begin ().sxDbgMsg([](auto _it){cout << "juhu:"<<_it).isValid ());
	//
	//   if (list.begin ().SX_DBG_MSG("juhu:"<<*iter).isValid ())
	//      cout << "valid\n";


	return 0;

	/*cout << list << endl;
	SX_CHECK (list.getSize() == 5, list);
	list.foreach ([](auto it) { cout << "*it=" << *it << endl; });


	SxList<int>::ConstIterator cit;
	SxList<int>::Iterator it;
	for (it = list.fromLast (); it != list.toFirst(); --it)  {
	cout << *it << endl;
	}


	list << 1 << 2 << 3 << SxRange<100,300,50>() << 1000 << 1001;
	cout << list << endl;
	cout << list.last () << endl;

	//   cout << SX_SEPARATOR;
	typedef SxList<int>::ConstIterator IntIterator;
	SxList<int> a = list.findAll ([](IntIterator it_)->bool { return *it_ < 250; });
	cout << SX_SEPARATOR;
	cout << a << endl;
	cout << SX_SEPARATOR;

	//   return 0;



	for (cit = list.begin(); cit != list.end(); ++cit) {
	cout << "it = " << *cit << endl;
	}

	list.foreach ([](SxList<int>::ConstIterator cit_)  {
	cout << *cit_ << endl;
	});

	// --- class iterator: {5, 4, 3, 2, 1, 0}
	for (long j : SxRange<5,0>()) { cout << "i=" << j << endl; }

	// --- for iterator loop
	typedef SxRange<13,-13> RangeType;
	RangeType range;
	RangeType::ConstIterator rangeIt;
	for (rangeIt = range.begin(); rangeIt != range.end(); ++rangeIt)  {
	cout << "rangeIt = " << *rangeIt << endl;
	}

	// --- sx::foreach
	sx::foreach (range.begin(), range.end(), [](RangeType::ConstIterator it_) {
	cout << "sx::foreach: " << *it_ << endl;
	});

	// --- member foreach
	range.foreach ([](RangeType::ConstIterator it_) {
	cout << "Range::foreach: " << *it_ << endl;
	});

	return 0;



	//   // --- foreach
	//   Range<0,100,10>().begin().foreach ([](Range<0,100,10>::ConstIterator it) {
	//      cout << *it << endl;
	//   });

	// SxList<int>::ConstIterator it;
	// for (it = i.begin(); it != i.end(); ++it)  {
	//    cout << *it << endl;
	// }


	return 0;  */
	/*

	SxList<int> ints;
	typedef SxList<int>::Iterator IntIt;
	ints << 1 << 2 << 3 << 4 << 5;

	cout << "sx::foreach" << endl;
	sx::foreach (ints.begin(), ints.end(), [](IntIt &it) {
	cout << *it << endl;
	});
	cout << SX_SEPARATOR;

	cout << "list.foreach" << endl;
	ints.foreach ([](IntIt it) { cout << *it << endl; });
	cout << SX_SEPARATOR;

	cout << "array.foreach" << endl;
	SxArray<int> intArray = ints;
	intArray.foreach ([](SxArray<int>::Iterator &it) { cout << *it << endl; });
	cout << SX_SEPARATOR;

	cout << "sx::find" << endl;
	if (sx::find (ints.begin(), ints.end(), 3))
	cout << "Found element" << endl;
	else
	cout << "Element not found" << endl;
	cout << SX_SEPARATOR;

	cout << "list.find" << endl;
	if (ints.find (3))  cout << "Found element" << endl;
	else                cout << "Element not found" << endl;
	cout << SX_SEPARATOR;

	cout << "sx::findCond" << endl;
	if (sx::findCond (ints.begin(), ints.end(),
	[](IntIt it)->bool {return *it > 2;}))
	cout << "Found element" << endl;
	else
	cout << "Element not found" << endl;
	cout << SX_SEPARATOR;

	cout << "list.findCond" << endl;
	if (ints.findCond ([](IntIt it)->bool { return *it > 2; }))
	cout << "Found element" << endl;
	else
	cout << "Element not found" << endl;
	cout << SX_SEPARATOR;

	cout << "Conditional results" << endl;
	SxList<int> res;
	ints.foreach ([&res](IntIt &it) { if (*it >= 2 && *it <= 4)  res << *it; });
	cout << res;
	cout << SX_SEPARATOR; */


	return 0;
}
