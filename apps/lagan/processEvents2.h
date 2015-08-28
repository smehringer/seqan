/*
 * processEvents.h
 *
 *  Created on: 22.06.2015
 *      Author: svenja
 */

#ifndef APPS_LAGAN_PROCESSEVENTS_H_
#define APPS_LAGAN_PROCESSEVENTS_H_

#include <seqan/stream.h>
#include <seqan/sequence.h>
#include <seqan/sequence_journaled.h>

using namespace seqan;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

enum DeltaEventType
{
	DELTA_EVENT_SNP = 0,
	DELTA_EVENT_DEL = 1,
	DELTA_EVENT_SV  = 2,
	DELTA_EVENT_INS = 3
};


//template<typename TValue>
struct DeltaEvent
{
	typedef Pair<unsigned, String<Dna5> > TPair;

	unsigned pos; // absolute position in reference sequence
	String<bool, Packed<> > seqs;
	DeltaEventType type;

	//values
	String<Dna5> ins;
	unsigned l_ins; // only for debugging reasons...!
	unsigned del;

	DeltaEvent(unsigned p, unsigned n, unsigned s, DeltaEventType t, String<Dna5> i, unsigned d):
		pos(p), type(t), ins(i), del(d)
	{
		resize(seqs, n, false);
		seqs[s] = true;
		l_ins = length(i);
	}

};

struct CompareByPosAndTypeLessThan_
{
	inline bool operator()(DeltaEvent const &lhs, DeltaEvent const &rhs)
	{
		if (lhs.pos == rhs.pos)
		{
			if (lhs.type == rhs.type)
			{
				if (lhs.ins == rhs.ins)
					return lhs.del < rhs.del;
				else
					return lhs.ins < rhs.ins;
			}
			return lhs.type < rhs.type;
		}
		return lhs.pos < rhs.pos;
	}
};

struct CompareByPosLessThan_
{
	inline bool operator()(DeltaEvent const &lhs, DeltaEvent const &rhs)
	{
		return lhs.pos < rhs.pos;
	}
};

struct DependentRegion
{
	typedef String<DeltaEvent> TDeltaEvents;
	typedef String<int> TScore;
	typedef StringSet<String<unsigned> > TDeps;

	TDeltaEvents records;
	TScore scores;
	TDeps dependencies;
	Pair<unsigned, int> bestScoringDeltaEvent;
};
// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================
int size(String<bool, Packed<> > & seqs)
{
	unsigned s = 0;
	for (unsigned i = 0; i < length(seqs); ++i)
		if (seqs[i] == true)
			++s;
	return s;
}

unsigned endPos(DeltaEvent & e)
{
	return e.pos + e.del; // deletion length
}

template < typename TString >
int deleteEntries(TString & string, String<unsigned> & dels)
{
	appendValue(dels, length(string)+1); // to avoid access error
	if (length(dels) != 0)
	{
		//assume dels are in ascending order
		TString tmp_string; // maybe resize before
		unsigned curr = 0;
		for (unsigned i = 0; i < length(string); ++i )
		{
			if (i != dels[curr])
				appendValue(tmp_string, string[i]);
			else
				++curr;
		}
		string = tmp_string;
	}
	return 0;
}

int scoreSNP(DependentRegion & dr, DeltaEvent & e, String<unsigned> & deps)
{
	// score = - size(e.seqs) + (length(e.seqs)-size(e.seqs))
	int score = - size(e.seqs);

	// snp:	nothing	 because another snp would not be affected by the change
	// del:		"	 because inside a deletion a snp-change doesn't matter
	// ins:		"	 because an insertion would just elongate by one character
	// sv :		"	 because of arguments from del and ins
	String<bool, Packed<> > tmp_seqs; // seqs that have at least one dependent delta event
	resize(tmp_seqs, length(e.seqs), false);
	tmp_seqs |= e.seqs;

	for (unsigned i = 0; i < length(deps); ++i)
		tmp_seqs |= dr.records[deps[i]].seqs; // so no sequence is subtracted twice

	score += (length(e.seqs)-size(tmp_seqs));
	return score;
}

int scoreDEL(DependentRegion & dr, DeltaEvent & e, String<unsigned> & deps)
{
	// subtract all deletions (one JE per deletion)
	int score = - size(e.seqs);

	// snp:(++score) because a snp would be included in the resulting insertion one JE must be added
	// del:(++score) because if another del is completely inside e the insertion causes one additional JE
	//                       else it forms a structural variant causing only one additional JE also
	// ins:          because an insertion just grows, no additional JE
	// sv :          because like an insertion the sv just alters without an additional JE

	String<bool, Packed<> > tmp_seqs; // seqs that have at least one dependent delta event
	resize(tmp_seqs, length(e.seqs), false);
	tmp_seqs |= e.seqs;

	// add one JE for each snp or deletion (except if deletion spans over e)
	for (unsigned i = 0; i < length(deps); ++i)
	{
		DeltaEvent & e2 = dr.records[deps[i]];
		if (e2.type == DELTA_EVENT_SNP) // snp
			score += size(e2.seqs);
		if (e2.type == DELTA_EVENT_DEL) // del
			if ((e.pos > e2.pos) | (endPos(e) < endPos(e2))) // if del_e2 does not include del_e
				score += size(e2.seqs);
		tmp_seqs |= e2.seqs;
		if (score > 0)
			return 0;
	}

	// add two JE's for every sequence that has no dependent delta event
	score += 2*(length(e.seqs)-size(tmp_seqs));

	return score;
}

int scoreINS(DependentRegion & dr, DeltaEvent & e, String<unsigned> & deps)
{
	// subtract all insertions (two JE per insertion)
	int score = - 2*size(e.seqs);

	// snp:(++score) because a snp would lead to a sv, one JE must be added
	// del:			 because del would just get bigger
	// ins:          because an insertion would lead to a sv (no additional JE)
	// sv :          because sv would just be altered
	String<bool, Packed<> > tmp_seqs; // seqs that have at least one dependent delta event
	resize(tmp_seqs, length(e.seqs), false);
	tmp_seqs |= e.seqs;

	// add one JE for each snp
	for (unsigned i = 0; i < length(deps); ++i)
	{
		DeltaEvent & e2 = dr.records[deps[i]];
		if (e2.type == DELTA_EVENT_SNP) // snp
			score += size(e2.seqs);
		tmp_seqs |= e2.seqs;
		if (score > 0)
			return 0;
	}

	// add one JE for every sequence that has no dependent delta event
	score += (length(e.seqs)-size(tmp_seqs));

	return score;
}

int scoreSV(DependentRegion & dr, DeltaEvent & e, String<unsigned> & deps)
{
	// subtract all sv's (two JE per sv)
	int score = - 2*size(e.seqs);
	//std::cout << e.del << " " << e.ins << "\n";
	// snp:(++score) because a snp would lead to a sv, one JE must be added
	// del:(++score) because del would lead to a sv (except if del fully includes event e)
	// ins:          because an insertion would lead to a sv (no additional JE)
	// sv :          because sv would just be altered
	String<bool, Packed<> > tmp_seqs; // seqs that have at least one dependent delta event
	resize(tmp_seqs, length(e.seqs), false);
	tmp_seqs |= e.seqs;

	// add one JE for each snp or delettion
	for (unsigned i = 0; i < length(deps); ++i)
	{
		DeltaEvent & e2 = dr.records[deps[i]];
		if ((e2.type == DELTA_EVENT_SNP)||(e2.type == DELTA_EVENT_DEL))// snp
			score += size(e2.seqs);
		tmp_seqs |= e2.seqs;
		if (score > 0)
			return 0;
	}

	// add two JE for every sequence that has no dependent delta event
	score += 2*(length(e.seqs)-size(tmp_seqs));

	if (score == 0) // if the change would not effect the number of journal entries, a smaller insertion is favored
		if (e.del < length(e.ins))
			--score;
	return score;
}

int computeEventScores(DependentRegion & dr)
{
	// score(DeltaEvent e) = if e would be merged into the reference seq
	// it is the number of journal entries that would emerge/diminish
	for (unsigned i = 0; i < length(dr.records); ++i)
	{
		DeltaEvent & e = dr.records[i];
		String<unsigned> & deps = dr.dependencies[i];
		int score;
		if (e.type == DELTA_EVENT_SNP)
			score = scoreSNP(dr, e, deps);
		else if (e.type == DELTA_EVENT_DEL)
			score = scoreDEL(dr, e, deps);
		else if (e.type == DELTA_EVENT_INS)
			score = scoreINS(dr, e, deps);
		else
			score = scoreSV(dr, e, deps);

		// minimum score is best -> lowest number of JE's
		if (score < dr.bestScoringDeltaEvent.i2)
		{
			Pair<unsigned, unsigned> best(i, score);
			dr.bestScoringDeltaEvent = best;
		}
	}

	return 0;
}


template<typename TSequence>
int mergeIntoRef(DeltaEvent & e, TSequence & ref)
{
	typedef typename Value<TSequence>::Type TSeqValue;
	typedef String<TSeqValue, Journaled<Alloc<>, SortedArray, Alloc<> > > TJournaledString;

	TJournaledString journal;
	setHost(journal, ref);

	if (e.type == DELTA_EVENT_SNP) // snp
		assignValue(journal, e.pos, e.ins[0]);
	else if (e.type == DELTA_EVENT_DEL) // del
		erase(journal, e.pos, endPos(e));
	else if (e.type == DELTA_EVENT_INS) // ins
		insert(journal, e.pos, e.ins);
	else // sv
	{
		erase(journal, e.pos, endPos(e));
		insert(journal, e.pos, e.ins);
	}

	flatten(journal);
	return 0;
}

template < typename TString, typename TInfix>
int getString(TString & s, TInfix & inf)
{
	typedef typename Iterator<TInfix>::Type TIterator;

	for(TIterator it = begin(inf); it != (end(inf)); goNext(it))
		appendValue(s, *it);
	return 0;
}

template<typename TRecords>
int applyOffset(TRecords & records, int offset, unsigned start)
{
	// assume records are sorted ascending by ref position
	for (unsigned i = start; i < length(records); ++i)
		records[i].pos += offset;
	return 0;
}

int updateDependencies(DependentRegion & dr)
{
	typedef Pair<unsigned, unsigned> TPair;

	clear(dr.dependencies);
	String<unsigned> dep_0 = "";
	appendValue(dr.dependencies, dep_0); // add string for first record

	String<TPair> dep;
	TPair pair(0, endPos(dr.records[0]));
	appendValue(dep, pair);

	for (unsigned i = 1; i < length(dr.records); ++i)
	{
		DeltaEvent & e = dr.records[i];
		String<unsigned> dep_i = "";

		String<unsigned> toDelete;
		for (unsigned d = 0; d < length(dep); ++d)
		{
			if (e.pos < dep[d].i2) // e and dep[d] are dependent of each other
			{
				appendValue(dr.dependencies[dep[d].i1], i);
				appendValue(dep_i, dep[d].i1);
			}
			else
			{
				appendValue(toDelete, d); // if current e.pos is not in range anymore it will never be again
			}
		}
		deleteEntries(dep, toDelete); //toDelete == 0?
		appendValue(dr.dependencies, dep_i);
		TPair next(i, endPos(e));
		appendValue(dep, next);
	}

	SEQAN_ASSERT(length(dr.dependencies)==length(dr.records));
	for (unsigned i = 0; i < length(dr.dependencies); ++i)
		for (unsigned j = 1; j < length(dr.dependencies[i]); ++j)
			SEQAN_ASSERT(dr.dependencies[i][j] > dr.dependencies[i][j-1]);

	return 0;
}

template<typename TList, typename TNeedle>
bool isIn(TList & list, TNeedle & needle)
{
	typedef typename Iterator<TList>::Type TIter;
	for (TIter it = begin(list); it != end(list); ++it)
	{
		if (*it == needle)
			return true;
	}
	return false;
}

// ------------ new
DeltaEventType determineType(unsigned del, String<Dna5> ins)
{
	DeltaEventType type;

	if (length(ins)==1 && del == 1)
		type = DELTA_EVENT_SNP;
	else if (length(ins)==0)
		type = DELTA_EVENT_DEL;
	else if (del == 0)
		type = DELTA_EVENT_INS;
	else
		type = DELTA_EVENT_SV;

	return type;
}

template<typename TSequence>
int processSingleEvent(DeltaEvent & d, DeltaEvent & e, TSequence const& ref)
{
	unsigned pos;
	DeltaEventType type;
	String<Dna5> ins = "";
	unsigned del = d.del;
	//int vp;

	if (d.pos <= e.pos)
	{
		pos = d.pos;
		append(ins, d.ins);
		SEQAN_ASSERT(ins == d.ins);
		del -= (endPos(d) - e.pos); // if to much is subtracted it will be corrected afterwards
	}
	else
	{
		pos = e.pos;
		String<Dna5> tmp_ins = infix(ref, e.pos, d.pos);
		append(ins, tmp_ins);
		append(ins, d.ins);
		del -= (endPos(d) - d.pos); // del = 0 always ?
	}

	if (endPos(d) <= endPos(e))
	{
		String<Dna5> tmp_ins = infix(ref, endPos(d), endPos(e));
		append(ins, tmp_ins);
	}
	else
	{
		del += (endPos(d) - endPos(e));
	}
	
	del += length(e.ins);

	//determine new type
	type = determineType(del, ins);

	// create new event and replace it with d
	DeltaEvent new_e = DeltaEvent(pos, length(e.seqs), 0, type, ins, del);
	new_e.seqs[0] = false; // revert initialization
	new_e.seqs |= d.seqs; // update seqs
	d = new_e; // assign event

	return 0;
}

DeltaEvent addEventToPrevious(DeltaEvent & e, DeltaEvent & p, int vp)
{
	unsigned pos = e.pos + vp;
	SEQAN_ASSERT(p.pos < pos);
	SEQAN_ASSERT((pos-p.pos) < length(p.ins)); // otherwise there would be no dependence...
	unsigned del = p.del;
	String<Dna5> ins = "";
	DeltaEventType type;

	// split p.ins into two at e.pos
	// insert e.ins there
	String<Dna5> ins1 = infix(p.ins, 0, (pos-p.pos));
	String<Dna5> ins2 = "";
	if ((pos-p.pos+e.del) < length(p.ins))
		ins2 = infix(p.ins, (pos-p.pos+e.del), length(p.ins));
	else
		del += (pos-p.pos+e.del) - length(p.ins);
	append(ins, ins1);
	append(ins, e.ins);
	append(ins, ins2);

	type = determineType(del, ins);

	DeltaEvent new_e = DeltaEvent(p.pos, length(p.seqs), 0, type, ins, del);
	new_e.seqs[0] = false; // revert initialization
	return new_e;
}

template<typename TSequence>
int processDR(DependentRegion & dr, TSequence & ref)
{
	computeEventScores(dr);
	bool better = dr.bestScoringDeltaEvent.i2 < 0;
	int offset = 0;

	while (better)
	{
		unsigned best_i = dr.bestScoringDeltaEvent.i1;
		DeltaEvent & best = dr.records[best_i];
		String<unsigned> dep_i = dr.dependencies[best_i];
		int tmp_offset = length(best.ins) - best.del;

		String<DeltaEvent> recs;
		String<int> vps; // store virtual positions in case of multiple events
		String<bool, Packed<> > unaffected_seqs;
		resize(unaffected_seqs, length(best.seqs), true);
		unaffected_seqs &= (~best.seqs);

		// inside a dependent region an event must not be dependant of every(!) other event
		// process dependent events
		for (unsigned d = 0; d < length(dep_i); ++d)
		{
			DeltaEvent & e = dr.records[dep_i[d]];
			unaffected_seqs &= (~e.seqs);

			// look for previous events on the same seq
			for (int r = length(recs)-1; r >= 0; --r)
			{
				String<bool, Packed<> > multiple = e.seqs & recs[r].seqs;
				if(!(testAllZeros(multiple)))
				{
					// event must be split up
					DeltaEvent new_e = addEventToPrevious(e, recs[r], vps[r]);
					new_e.seqs = multiple;
					e.seqs &= (~multiple); // remove multiple
					recs[r].seqs &= (~multiple); // remove multiple
					if (testAllZeros(recs[r].seqs))
					{
						assignValue(recs, r, new_e); // can be replaced
						vps[r]+= ((length(e.ins)-e.del));
					}
					else
					{
						appendValue(recs, new_e);
						appendValue(vps, vps[r]+(length(e.ins)-e.del));
					}
				}
				if (testAllZeros(e.seqs))
					break;
			}
			if (!(testAllZeros(e.seqs)))
			{
				appendValue(vps, (length(e.ins)-e.del));
				processSingleEvent(e, best, ref);
				appendValue(recs, e);
			}
		}

		// add independent records to recs
		// apply offset to independent events after best
		for (unsigned i = 0; i < length(dr.records); ++i)
		{
			if (!(isIn(dep_i, i)) && (i!=best_i))
			{
				DeltaEvent & ev = dr.records[i];
				if (ev.pos >= endPos(best))
					ev.pos += tmp_offset;
				appendValue(recs, ev);
			}
		}

		// now create new event (opposite of best) for all unaffected sequences
		unsigned del = length(best.ins);
		String<Dna5> ins = infix(ref, best.pos, endPos(best));
		DeltaEventType type = determineType(del, ins);
		DeltaEvent new_e = DeltaEvent(best.pos, length(best.seqs), 0, type, ins, del);
		new_e.seqs = unaffected_seqs;
		append(recs, new_e);

		mergeIntoRef(best, ref);

		//update dr and look for other events to merged into ref
		dr.records = recs;
		dr.bestScoringDeltaEvent.i1 = 0;
		dr.bestScoringDeltaEvent.i2 = 0;

		sort(dr.records, CompareByPosLessThan_());
		updateDependencies(dr);
		computeEventScores(dr);
		better = dr.bestScoringDeltaEvent.i2 < 0;
		offset += tmp_offset;
	}

	return offset;
}



bool isEqual(DeltaEvent & lhs, DeltaEvent & rhs)
{
	if (lhs.pos == rhs.pos)
		if (lhs.type == rhs.type)
			if (lhs.del == rhs.del)
				return lhs.ins == rhs.ins;
	return false;
}

template < typename TPair>
int updateDR(DependentRegion & dr, DeltaEvent & e, TPair & dep)
{
	unsigned record_num = length(dr.records); // index in dependent region

	// check which events are influenced by delta event e
	String<unsigned> toDelete = "";
	String<unsigned> toAdd = "";
	for (unsigned d = 0; d < length(dep); ++d)
	{
		if (e.pos < dep[d].i2) // e and dep[d] are dependent of each other
		{
			appendValue(dr.dependencies[dep[d].i1], record_num);
			appendValue(toAdd, dep[d].i1);
		}
		else
		{
			appendValue(toDelete, d); // if current e.pos is not in range anymore it will never be again
		}
	}
	deleteEntries(dep, toDelete); //toDelete == 0?

	appendValue(dr.records, e);
	appendValue(dr.dependencies, toAdd);
	return 0;
}

template<typename TDeltaEvents>
unsigned getNextDR(DependentRegion & dr, TDeltaEvents & records, unsigned start)
{
	SEQAN_ASSERT(start < length(records));
	typedef Pair<unsigned, unsigned> TPair;
	typedef String<bool, Packed<> > TPacked;

	TPacked allZeros;
	resize(allZeros, length(records[0].seqs), false);
	String<TPair> dep;

	//DeltaEvent & first = records[start];
	// initialize first dependent region
	appendValue(dr.records, records[start]);
	String<unsigned> d = ""; // no dependencies for first entry yet
	appendValue(dr.dependencies, d);

	unsigned drb = endPos(records[start]); // dependent Region Border
	TPair pair(0, drb);
	appendValue(dep, pair);


	// parse records from left to right
	// look for merging possibilities and define dependent regions
	for (unsigned i = start+1; i < length(records); ++i)
	{
		DeltaEvent & e = records[i];

		if (isEqual(e, records[i-1]))
		{
			dr.records[length(dr.records)-1].seqs |= e.seqs;
			e.seqs &= allZeros;
		}
		else if (e.pos < drb) // e is inside dependent region dr
		{
			unsigned record_num = length(dr.records); // index in dependent region

			updateDR(dr, e, dep);

			unsigned b = endPos(e);
			if (b > drb)
				drb = b;
			TPair pair(record_num, b);
			appendValue(dep, pair);
		}
		else // end of dependent region
		{
			return (i); // next region starts with i
		}
	}

	return length(records);
}

// ----------------------------------------------------------------------------
// Function processDeltaEventsOnReference()
// ----------------------------------------------------------------------------
//main
template<typename TDeltaEvents, typename TSequence>
int processDeltaEventsOnReference(TDeltaEvents & records, TSequence & ref)
{
	String<DeltaEvent> recs;
	sort(records, CompareByPosAndTypeLessThan_()); // sort by reference position(first) and type(second) and value(third)

	unsigned start = 0;
	unsigned end = length(records);
	while(start < end)
	{
		DependentRegion dr;
		unsigned new_start = getNextDR(dr, records, start);
		int offset = processDR(dr, ref); // inside here, no record is deleted

		// copy records from dr to new list because number of records varys now
		for(unsigned i = 0; i < length(dr.records); ++i)
			appendValue(recs, dr.records[i]);

		if(offset!=0)
			applyOffset(records, offset, new_start);
		start = new_start;
		end = length(records); // update because number of records might change!
	}

	records = recs;

	return 0;
}
#endif /* APPS_LAGAN_PROCESSEVENTS_H_ */
