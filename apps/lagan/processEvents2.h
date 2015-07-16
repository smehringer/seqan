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
//template<typename TValue>
struct DeltaEvent
{
	typedef Pair<unsigned, String<Dna5> > TPair;

	unsigned pos; // absolute position in reference sequence
	String<bool, Packed<> > seqs;
	unsigned type; // 0 = SNP, Del = 1, Ins = 2, SV = 3

	//values
	String<Dna5> ins;
	unsigned del;

	DeltaEvent(unsigned p, unsigned s, unsigned n, unsigned t, String<Dna5> i, unsigned d) :
		pos(p), type(t), ins(i), del(d),
	{
		resize(seqs, n, false);
		seqs[s] = true;
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
				if ((lhs.type == 0) || (lhs.type == 2))
				{
					return lhs.ins < rhs.ins;
				}
				else if (lhs.type == 1)
				{
					return lhs.del < rhs.del;
				}
				else // lhs.type = 3
				{
					if (lhs.ins == rhs.ins)
						return lhs.del < rhs.del;
					else
						return lhs.ins < rhs.ins;
				}
			}
			else
				return lhs.type < rhs.type;
		}
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

int combineOR(String<bool, Packed<> > & seqs1, String<bool, Packed<> > & seqs2)
{
	SEQAN_ASSERT_EQ(length(seqs1), length(seqs2));
	for (unsigned i = 0; i < length(seqs1); ++i)
	{
		bool b = seqs1[i] or seqs2[i];
		seqs1[i] = b;
	}
	return 0;
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
	int score = (int)length(e.seqs) - 2*size(e.seqs);

	// snp:(--score) because another snp would not be affected by the change
	// del:(--score) because inside a deletion a snp-change doesn't matter
	// ins:(--score) because an insertion would just elongate by one character
	// sv :(--score) because of arguments from del and ins
	// -> substract all those unaffected sequences
	String<bool, Packed<> > tmp_seqs;
	resize(tmp_seqs, length(e.seqs));
	for(unsigned i = 0; i < length(e.seqs); ++i)
		tmp_seqs[i] = false;

	for (unsigned i = 0; i < length(deps); ++i)
		combineOR(tmp_seqs, dr.records[deps[i]].seqs); // so no sequence is subtracted twice

	score -= size(tmp_seqs);
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
	resize(tmp_seqs, length(e.seqs));
	for(unsigned i = 0; i < length(e.seqs); ++i)
		tmp_seqs[i] = false;
	combineOR(tmp_seqs, e.seqs);

	// add one JE for each snp or deletion (except if deletion spans over e)
	for (unsigned i = 0; i < length(deps); ++i)
	{
		DeltaEvent & e2 = dr.records[deps[i]];
		if (e2.type == 0) // snp
			score += size(e2.seqs);
		if (e2.type == 1) // del
			if ((e.pos > e2.pos) | (endPos(e) < endPos(e2))) // if del_e2 does not include del_e
				score += size(e2.seqs);
		combineOR(tmp_seqs, e2.seqs);
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
	resize(tmp_seqs, length(e.seqs));
	for(unsigned i = 0; i < length(e.seqs); ++i)
		tmp_seqs[i] = false;
	combineOR(tmp_seqs, e.seqs);

	// add one JE for each snp
	for (unsigned i = 0; i < length(deps); ++i)
	{
		DeltaEvent & e2 = dr.records[deps[i]];
		if (e2.type == 0) // snp
			score += size(e2.seqs);
		combineOR(tmp_seqs, e2.seqs);
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

	// snp:(++score) because a snp would lead to a sv, one JE must be added
	// del:(++score) because del would lead to a sv (except if del fully includes event e)
	// ins:          because an insertion would lead to a sv (no additional JE)
	// sv :          because sv would just be altered
	String<bool, Packed<> > tmp_seqs; // seqs that have at least one dependent delta event
	resize(tmp_seqs, length(e.seqs));
	for(unsigned i = 0; i < length(e.seqs); ++i)
		tmp_seqs[i] = false;
	combineOR(tmp_seqs, e.seqs);

	// add one JE for each snp or delettion
	for (unsigned i = 0; i < length(deps); ++i)
	{
		DeltaEvent & e2 = dr.records[deps[i]];
		if (e2.type == 0) // snp
			score += size(e2.seqs);
		if (e2.type == 1) // del
			if ((e.pos > e2.pos) | (endPos(e) < endPos(e2))) // if del_e2 does not include del_e
				score += size(e2.seqs);
		combineOR(tmp_seqs, e2.seqs);
		if (score > 0)
			return 0;
	}

	// add two JE for every sequence that has no dependent delta event
	score += 2*(length(e.seqs)-size(tmp_seqs));

	if (score == 0) // if the change would not effect the number of journal entries, a smaller insertion is favoured
		if (length(e2.ins) < length(e.ins))
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
		if (e.type == 0)
			score = scoreSNP(dr, e, deps);
		else if (e.type == 1)
			score = scoreDEL(dr, e, deps);
		else if (e.type == 2)
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

	if (e.type == 0) // snp
		assignValue(journal, e.pos, e.ins[0]);
	else if (e.type == 1) // del
		erase(journal, e.pos, endPos(e));
	else if (e.type == 2) // ins
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

	SEQAN_ASSERT_EQ(length(dr.records), length(dr.dependencies));
	for (unsigned i = 1; i < length(dr.records); ++i) // maybe on the outside
		clear(dr.dependencies[i]);

	String<TPair> dep;
	TPair pair(0, endPos(dr.records[0]));
	appendValue(dep, pair);

	for (unsigned i = 1; i < length(dr.records); ++i)
	{
		DeltaEvent & e = dr.records[i];

		String<unsigned> toDelete;
		String<unsigned> toAdd;
		for (unsigned d = 0; d < length(dep); ++d)
		{
			if (e.pos < dep[d].i2) // e and dep[d] are dependent of each other
			{
				appendValue(dr.dependencies[dep[d].i1], i);
				appendValue(toAdd, dep[d].i1);
			}
			else
			{
				appendValue(toDelete, d); // if current e.pos is not in range anymore it will never be again
			}
		}
		deleteEntries(dep, toDelete); //toDelete == 0?
		dr.dependencies[i] = toAdd;
		TPair next(i, endPos(e));
		appendValue(dep, next);
	}

	return 0;
}

template<typename TSequence>
int updateSNP(DependentRegion & dr, DeltaEvent & e, unsigned i, TSequence & ref)
{
//	std::cout << "included SNP at " << e.pos << " into ref "<< std::endl;
	typedef Pair<unsigned, String<Dna5> > Tsv;

	String<unsigned> & deps = dr.dependencies[i];

	String<bool, Packed<> > tmp_seqs; // seqs that have at least one dependent delta event
	resize(tmp_seqs, length(e.seqs));
	for(unsigned j = 0; j < length(e.seqs); ++j)
		tmp_seqs[j] = false;
	combineOR(tmp_seqs, e.seqs);

	// alter sequences with dependent events
	// only insertion needs an update
	for (unsigned d = 0; d < length(deps); ++d)
	{
		DeltaEvent & dep_e = dr.records[deps[d]];
		if (dep_e.type == 2) // convert ins to sv
		{
			dep_e.type = 3;
			Tsv sv(1, dep_e.ins);
			dep_e.sv = sv;
			dep_e.ins = 0;
		}
		combineOR(tmp_seqs, dep_e.seqs);
	}

	// create new snp for every sequence without an dependent event
	DeltaEvent new_e = DeltaEvent(e.pos, 0, length(e.seqs), ref[e.pos]);
	for (unsigned s = 0; s < length(e.seqs); ++s)
		new_e.seqs[s] = not tmp_seqs[s];

	//erase(dr.records, i)
	//std::cout << "Before Assignment: "<< dr.records[i].pos << " " << dr.records[i].snp << " seqs: "<< dr.records[i].seqs[0]<< dr.records[i].seqs[1]<< dr.records[i].seqs[2]<< dr.records[i].seqs[3] << std::endl;
	assignValue(dr.records, i, new_e);
	//std::cout << "After Assignment: "<< dr.records[i].pos << " " << dr.records[i].snp << " seqs: "<< dr.records[i].seqs[0]<< dr.records[i].seqs[1]<< dr.records[i].seqs[2]<< dr.records[i].seqs[3] << std::endl;
	//std::cout << "Should be.......: "<< dr.records[i].pos << " " << ref[e.pos] <<  std::endl;
	return 0;
}

template<typename TSequence>
int updateDEL(DependentRegion & dr, DeltaEvent & e, unsigned i, TSequence & ref)
{
	std::cout << "included Deletion at " << e.pos << " into ref "<< std::endl;
	typedef Pair<unsigned, String<Dna5> > Tsv;

	String<unsigned> & deps = dr.dependencies[i];

	String<bool, Packed<> > tmp_seqs; // seqs that have at least one dependent delta event
	resize(tmp_seqs, length(e.seqs));
	for(unsigned j = 0; j < length(e.seqs); ++j)
		tmp_seqs[j] = false;
	combineOR(tmp_seqs, e.seqs);

	// alter sequences with dependent events
	for (unsigned d = 0; d < length(deps); ++d)
	{
		DeltaEvent & dep_e = dr.records[deps[d]];
		if (dep_e.type == 0) // snp -> insertion
		{
			String<Dna5> ins;
			Infix<Dna5String>::Type infix(ref, e.pos, getBorder(e));
			getString(ins, infix);
			assignValue(ins, abs(dep_e.pos-e.pos), dep_e.snp);
			dep_e.type = 2;
			dep_e.ins = ins;
			dep_e.snp = 0;
			dep_e.pos = e.pos;
		}
		else if (dep_e.type == 1) // del
		{
			if ((dep_e.pos >= e.pos) & (getBorder(dep_e) <= getBorder(e)))
			{
				String<Dna5> ins;
				Infix<Dna5String>::Type infix1(ref, e.pos, dep_e.pos);
				Infix<Dna5String>::Type infix2(ref, getBorder(dep_e), getBorder(e));
				getString(ins, infix1);
				getString(ins, infix2);
				dep_e.type = 2;
				dep_e.ins = ins;
				dep_e.del = 0;
				dep_e.pos = e.pos;
			}
			else if ((dep_e.pos >= e.pos) & (getBorder(dep_e) > getBorder(e)))
			{
				String<Dna5> ins;
				Infix<Dna5String>::Type infix(ref, e.pos, dep_e.pos);
				getString(ins, infix);
				unsigned del = getBorder(dep_e) - getBorder(e);
				dep_e.type = 3;
				Tsv sv(del, ins);
				dep_e.sv = sv;
				dep_e.del = 0;
				dep_e.pos = e.pos;
			}
			else if ((dep_e.pos < e.pos) & (getBorder(dep_e) <= getBorder(e)))
			{
				String<Dna5> ins;
				Infix<Dna5String>::Type infix(ref, getBorder(dep_e), getBorder(e));
				getString(ins, infix);
				unsigned del = e.pos - dep_e.pos;
				dep_e.type = 3;
				Tsv sv(del, ins);
				dep_e.sv = sv;
				dep_e.del = 0;
			}
			else
			{
				dep_e.del -= e.del;
			}
		}
		else if (dep_e.type == 2) // insertion grows
		{
			String<Dna5> ins;
			Infix<Dna5String>::Type infix1(ref, e.pos, dep_e.pos);
			Infix<Dna5String>::Type infix2(ref, dep_e.pos, getBorder(e));
			getString(ins, infix1);
			append(ins, dep_e.ins);
			getString(ins, infix2);
			dep_e.ins = ins;
		}
		else if (dep_e.type == 3) // alter sv
		{
			if ((dep_e.pos >= e.pos) & (getBorder(dep_e) <= getBorder(e))) // sv -> insertion
			{
				String<Dna5> ins;
				Infix<Dna5String>::Type infix1(ref, e.pos, dep_e.pos);
				Infix<Dna5String>::Type infix2(ref, getBorder(dep_e), getBorder(e));
				getString(ins, infix1);
				append(ins, dep_e.sv.i2);
				getString(ins, infix2);
				dep_e.type = 2;
				dep_e.ins = ins;
				dep_e.pos = e.pos;
				//dep_e.sv = 0; todo:: how to set no Null ?
			}
			else if ((dep_e.pos >= e.pos) & (getBorder(dep_e) > getBorder(e)))
			{
				String<Dna5> ins;
				Infix<Dna5String>::Type infix(ref, e.pos, dep_e.pos);
				getString(ins, infix);
				append(ins, dep_e.sv.i2);
				unsigned del = getBorder(dep_e) - getBorder(e);
				Tsv sv(del, ins);
				dep_e.sv = sv;
				dep_e.pos = e.pos;
			}
			else if ((dep_e.pos < e.pos) & (getBorder(dep_e) <= getBorder(e)))
			{
				String<Dna5> ins = dep_e.sv.i2;
				Infix<Dna5String>::Type infix(ref, getBorder(dep_e), getBorder(e));
				getString(ins, infix);
				unsigned del = e.pos - dep_e.pos;
				Tsv sv(del, ins);
				dep_e.sv = sv;
			}
			else
			{
				dep_e.sv.i1 -= e.del;
			}
			// else deletion of dep_e covers all of e so no change occurs
		}
		combineOR(tmp_seqs, dep_e.seqs);
	}

	// create new insertion for every sequence without an dependent event
	String<Dna5> ins;
	Infix<Dna5String>::Type infix(ref, e.pos, getBorder(e));
	getString(ins, infix);
	DeltaEvent new_e = DeltaEvent(e.pos, 0, length(e.seqs), ins);
	for (unsigned s = 0; s < length(e.seqs); ++s)
		new_e.seqs[s] = not tmp_seqs[s];

	//erase(dr.records, i);
	assignValue(dr.records, i, new_e);

	return 0;
}

int updateINS(DependentRegion & dr, DeltaEvent & e, unsigned i)
{
	std::cout << "included Insertion at " << e.pos << " into ref "<< std::endl;
	typedef Pair<unsigned, String<Dna5> > Tsv;

	String<unsigned> & deps = dr.dependencies[i];

	String<bool, Packed<> > tmp_seqs; // seqs that have at least one dependent delta event
	resize(tmp_seqs, length(e.seqs));
	for(unsigned j = 0; j < length(e.seqs); ++j)
		tmp_seqs[j] = false;
	combineOR(tmp_seqs, e.seqs);

	// alter sequences with dependent events
	for (unsigned d = 0; d < length(deps); ++d)
	{
		DeltaEvent & dep_e = dr.records[deps[d]];

		if (dep_e.type == 0) // snp -> sv
		{
			String<Dna5> ins;
			appendValue(ins, dep_e.snp);
			Tsv sv(length(e.ins)+1, ins);
			dep_e.type = 3;
			dep_e.sv = sv;
			dep_e.snp = 0;
			dep_e.pos = e.pos;
		}
		else if (dep_e.type == 1) // del : if insertion lies within a deletion, the deletion just gets bigger
		{
			dep_e.del += length(e.ins);
		}
		else if (dep_e.type == 2) // insertion -> sv
		{
			Tsv sv(length(e.ins), dep_e.ins);
			dep_e.type = 3;
			dep_e.sv = sv;
			dep_e.ins = 0;
		}
		else if (dep_e.type == 3) // sv : if insertion lies within a deletion of a sv, the deletion just gets bigger
		{
			dep_e.sv.i1 += length(e.ins);
		}
		combineOR(tmp_seqs, dep_e.seqs);
	}

	// create new deletion for every sequence without an dependent event
	unsigned l = length(e.ins);
	DeltaEvent new_e = DeltaEvent(e.pos, 0, length(e.seqs), l);
	for (unsigned s = 0; s < length(e.seqs); ++s)
		new_e.seqs[s] = not tmp_seqs[s];

	//erase(dr.records, i);
	assignValue(dr.records, i, new_e);

	return 0;
}

template<typename TSequence>
int updateSV(DependentRegion & dr, DeltaEvent & e, unsigned i, TSequence & ref)
{
	std::cout << "included structural variant at " << e.pos << " into ref "<< std::endl;
	typedef Pair<unsigned, String<Dna5> > Tsv;

	String<unsigned> & deps = dr.dependencies[i];

	String<bool, Packed<> > tmp_seqs; // seqs that have at least one dependent delta event
	resize(tmp_seqs, length(e.seqs));
	for(unsigned j = 0; j < length(e.seqs); ++j)
		tmp_seqs[j] = false;
	combineOR(tmp_seqs, e.seqs);

	// alter sequences with dependent events
	for (unsigned d = 0; d < length(deps); ++d)
	{
		DeltaEvent & dep_e = dr.records[deps[d]];

		if (dep_e.type == 0) // snp -> sv
		{
			String<Dna5> ins;
			Infix<Dna5String>::Type infix(ref, e.pos, getBorder(e));
			getString(ins, infix);
			assignValue(ins, abs(dep_e.pos-e.pos), dep_e.snp);
			unsigned del = length(e.sv.i2);
			dep_e.type = 3;
			Tsv sv(del, ins);
			dep_e.sv = sv;
			dep_e.snp = 0;
			dep_e.pos = e.pos;
		}
		else if (dep_e.type == 1) // del
		{
			if ((dep_e.pos >= e.pos) & (getBorder(dep_e) <= getBorder(e))) // del -> sv
			{
				String<Dna5> ins;
				Infix<Dna5String>::Type infix1(ref, e.pos, dep_e.pos);
				Infix<Dna5String>::Type infix2(ref, getBorder(dep_e), getBorder(e));
				getString(ins, infix1);
				getString(ins, infix2);
				dep_e.type = 3;
				Tsv sv(e.sv.i1, ins);
				dep_e.sv = sv;
				dep_e.del = 0;
				dep_e.pos = e.pos;
			}
			else if ((dep_e.pos >= e.pos) & (getBorder(dep_e) > getBorder(e))) // del -> sv
			{
				String<Dna5> ins;
				Infix<Dna5String>::Type infix(ref, e.pos, dep_e.pos);
				getString(ins, infix);
				unsigned del = getBorder(dep_e) - getBorder(e);
				dep_e.type = 3;
				Tsv sv(del, ins);
				dep_e.sv = sv;
				dep_e.del = 0;
				dep_e.pos = e.pos;
			}
			else if ((dep_e.pos < e.pos) & (getBorder(dep_e) <= getBorder(e))) // del -> sv
			{
				String<Dna5> ins;
				Infix<Dna5String>::Type infix(ref, getBorder(dep_e), getBorder(e));
				getString(ins, infix);
				unsigned del = e.pos - dep_e.pos + length(e.sv.i2);
				dep_e.type = 3;
				Tsv sv(del, ins);
				dep_e.sv = sv;
				dep_e.del = 0;
			}
			else // del fully includes sv
			{
				dep_e.del += ( length(e.sv.i2) - e.sv.i1);
			}
		}
		else if (dep_e.type == 2) // ins -> sv
		{
			String<Dna5> ins;
			Infix<Dna5String>::Type infix1(ref, e.pos, dep_e.pos);
			Infix<Dna5String>::Type infix2(ref, dep_e.pos, getBorder(e));
			getString(ins, infix1);
			append(ins, dep_e.ins);
			getString(ins, infix2);
			std::cout << ins << std::endl;
			Tsv sv(length(e.sv.i2), ins);
			dep_e.type = 3;
			dep_e.sv = sv;
			dep_e.ins = 0;
			dep_e.pos = e.pos;
		}
		else if (dep_e.type == 3) // alter sv
		{
			if ((dep_e.pos >= e.pos) & (getBorder(dep_e) <= getBorder(e)))
			{
				String<Dna5> ins;
				Infix<Dna5String>::Type infix1(ref, e.pos, dep_e.pos);
				Infix<Dna5String>::Type infix2(ref, getBorder(dep_e), getBorder(e));
				getString(ins, infix1);
				append(ins, dep_e.sv.i2);
				getString(ins, infix2);
				dep_e.sv.i1 = length(e.sv.i2);
				dep_e.sv.i2 = ins;
				dep_e.pos = e.pos;
			}
			else if ((dep_e.pos >= e.pos) & (getBorder(dep_e) > getBorder(e)))
			{
				String<Dna5> ins;
				Infix<Dna5String>::Type infix(ref, e.pos, dep_e.pos);
				getString(ins, infix);
				append(ins, dep_e.sv.i2);
				unsigned del = length(e.sv.i2) + getBorder(dep_e) - getBorder(e);
				dep_e.sv.i2 = ins;
				dep_e.sv.i1 = del;
				dep_e.pos = e.pos;
			}
			else if ((dep_e.pos < e.pos) & (getBorder(dep_e) <= getBorder(e)))
			{
				String<Dna5> ins = dep_e.sv.i2;
				Infix<Dna5String>::Type infix(ref, getBorder(dep_e), getBorder(e));
				getString(ins, infix);
				unsigned del = length(e.sv.i2) + e.pos - dep_e.pos;
				dep_e.sv.i1 = del;
				dep_e.sv.i2 = ins;
			}
			else
			{
				dep_e.sv.i1 += (length(e.sv.i2) - e.sv.i1);
			}
		}
		combineOR(tmp_seqs, dep_e.seqs);
	}

	// create new sv for every sequence without an dependent event
	String<Dna5> ins;
	Infix<Dna5String>::Type infix(ref, e.pos, getBorder(e));
	getString(ins, infix);
	Tsv sv(length(e.sv.i2), ins);
	DeltaEvent new_e = DeltaEvent(e.pos, 0, length(e.seqs), sv);
	for (unsigned s = 0; s < length(e.seqs); ++s)
		new_e.seqs[s] = not tmp_seqs[s];

	//erase(dr.records, i);
	assignValue(dr.records, i, new_e);

	return 0;
}

// ----------------------------------------------------------------------------
// Function processDR()
// ----------------------------------------------------------------------------
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
		DeltaEvent best_copy = dr.records[best_i];
		int tmp_offset = 0;

		// update other sequences (delta map)
		if (best.type == 0) // snp
		{
			updateSNP(dr, best, best_i, ref);
		}
		else if (best.type == 1) // del
		{
			tmp_offset = - best.del;
			updateDEL(dr, best, best_i, ref);
		}
		else if (best.type == 2) // ins
		{
			tmp_offset = length(best.ins);
			updateINS(dr, best, best_i);
		}
		else // sv
		{
			tmp_offset = length(best.sv.i2) - best.sv.i1;
			updateSV(dr, best, best_i, ref);
		}

		// include best event into reference
		mergeIntoRef(best_copy, ref);

		sort(dr.records, CompareByPosAndTypeLessThan_());

		// update tmp_offset on unaffected delta events in dr
		String<unsigned> dep_i = dr.dependencies[best_i];
		unsigned start;
		if (length(dep_i)== 0)
			start = best_i+1;
		else
		{
			unsigned last_dep = dep_i[length(dep_i)-1];
			if (last_dep < best_i)
				start = best_i +1 ;
			else
				start = last_dep +1;
		}
		applyOffset(dr.records, tmp_offset, start);

		// update dependencies
		updateDependencies(dr);

		// compute scores again
		dr.bestScoringDeltaEvent.i2 = 0;

		computeEventScores(dr);

		// set better again
		better = dr.bestScoringDeltaEvent.i2 < 0;
		offset += tmp_offset;
	}
	return offset;
}




// ------------ new
template<typename TSeq>
String<Dna5> myInfix(TSeq & seq, unsigned s, unsigned e)
{
	// makes a copy of seq from s (inklusive) to e (exklusive)

	String<Dna5> inf;
	if(s >= e)
		return inf; // todo:: return empty string?
	for (unsigned i = s; i < e; ++i)
		appendValue(inf, seq[i]);
	return inf;
}

int renewSeqs(DeltaEvent & e, unsigned s)
{
	for (unsigned i = 0; i < length(e.seqs); ++i)
		e.seqs = false;
	e.seqs[s] == true;
}

template<typename TSequence>
int processSingleEvent(DeltaEvent & d, DeltaEvent const& e, TSequence const& ref)
{
	unsigned pos;
	unsigned type;
	String<Dna5> ins = "";
	unsigned del = d.del;
	if (d.pos <= e.pos)
	{
		pos = d.pos;
		append(ins, d.ins);
		SEQAN_ASSERT(ins, d.ins);
		del -= (endPos(d) - e.pos); // if to much is subtracted it will be corrected afterwards
	}
	elseString<bool, Packed<> >
	{
		pos = e.pos;
		tmp_ins = myInfix(ref, e.pos, d.pos);
		append(ins, tmp_ins); 
		del -= (endPos(d) - d.pos); // del = 0 always ?
	}

	if (endPos(d) <= endPos(e))
	{
		tmp_ins = myInfix(ref, endPos(d), endPos(e);
		apennd(ins, tmp_ins);
	}
	else
	{
		del += (endPos(d) - endPos(e));
	}
	
	del += length(e.ins);

	//determine new type
	if (length(ins)==1 && del == 1)
		type = 0;
	else if (length(ins)==0)
		type = 1;
	else if (del == 0)
		type = 2;
	else
		type = 3;

	// create new event and replace it with d
	DeltaEvent new_e = DeltaEvent(pos, length(e.seqs), 0, type, ins, del);
	new_e.seqs[0] = false; // revert initialization
	new_e |= d.seqs; // update seqs
	d = new_e; // assign event

	return 0;
}

template<typename TSequence>
int processDRnew(DependentRegion & dr, TSequence & ref)
{
        computeEventScores(dr);
        bool better = dr.bestScoringDeltaEvent.i2 < 0;
        int offset = 0;
        while (better)
	{
		unsigned best_i = dr.bestScoringDeltaEvent.i1;
                DeltaEvent & best = dr.records[best_i];
                DeltaEvent best_copy = dr.records[best_i];
		String<unsigned> dep_i = dr.dependencies[best_i];
                int tmp_offset = 0;
		
		String<bool, Packed<> > unaffected_seqs;
		resize(unaffected_seqs, length(best.seqs), true);
		String<DeltaEvents> recs;

		// insige a dependent region an event must not be dependant of every(!) other event
		// process dependent events
		for (unsigned d = 0; d < length(dep_i); ++d)
		{
			DeltaEvent & e = records[dep_i[d]];
			// verunde unaffected_seqs and negative of e.seqs

			// look for previous events on the same seq
			for (int r = length(recs)-1; r >= 0; --r)
			{
				String<bool, Packed<> > multiple = e.seqs & recs[r];
				if(!(testAllZeros(multiple)))
				{
					// event must be split up
					new_e = addEventToPrevious(e, dr.records[dep_i[r]]);
					new_e.seqs = multiple;
					// e.seqs &= !(multiple); // remove multiple
					// recs[r].seqs &= !(multiple); // remove multiple
					if (testAllZeros(recs[r].seqs))
						assignValue(recs, r, new_e); // can be replaced
					else
						appendValue(recs, new_e);
				}
				if (testAllZeros(e.seqs))
					break;
			}
			if (!(testAllZeros(e.seqs)))
			{
				processSingleEvent(e, best, ref);
				appendValue(recs, e);
			}
		}
		
		// apply offset to independant events after best
	}

	return offset;

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

bool isEqual(DeltaEvent & lhs, DeltaEvent & rhs)
{
	if (lhs.pos == rhs.pos)
	{
		if (lhs.type == rhs.type)
		{
			if (lhs.type == 0)
				return lhs.snp == rhs.snp;
			else if (lhs.type == 1)
				return lhs.del == rhs.del;
			else if (lhs.type == 2)
				return lhs.ins == rhs.ins;
			else // lhs.type = 3
				return lhs.sv == rhs.sv; // pair ?
		}
	}
	return false;
}

template < typename TPair>
int updateDR(DependentRegion & dr, DeltaEvent & e, TPair & dep)
{
	unsigned record_num = length(dr.records); // index in dependent region

	// check which events are influenced by delta event e
	String<unsigned> toDelete;
	String<unsigned> toAdd;
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

int clear(DependentRegion & dr)
{
	clear(dr.records);
	clear(dr.scores);
	clear(dr.dependencies);
	Pair<unsigned, unsigned> pair(0,0);
	dr.bestScoringDeltaEvent = pair;

	return 0;
}

template<typename TDeltaEvents>
unsigned getNextDR(DependentRegion & dr, String<unsigned> & recs_to_delete, TDeltaEvents & records, unsigned start)
{
	SEQAN_ASSERT(start < length(records));
	typedef Pair<unsigned, unsigned> TPair;

	String<TPair> dep;

	//DeltaEvent & first = records[start];
	// initialize first dependent region
	appendValue(dr.records, records[start]);
	String<unsigned> d; // no dependencies for first entry yet
	appendValue(dr.dependencies, d);

	unsigned drb = getBorder(records[start]); // dependent Region Border
	TPair pair(0, drb);
	appendValue(dep, pair);


	// parse records from left to right
	// look for merging possibilities and define dependent regions
	//String<unsigned> recs_to_delete;
	for (unsigned i = start+1; i < length(records); ++i)
	{
		DeltaEvent & e = records[i];

		if (isEqual(e, records[i-1]))
		{
			// verodere sequence information und loesche e aus records
			combineOR(dr.records[length(dr.records)-1].seqs, e.seqs);
			appendValue(recs_to_delete, i);
		}

		else if (e.pos < drb) // e is inside dependent region dr
		{
			unsigned record_num = length(dr.records); // index in dependent region

			updateDR(dr, e, dep);

			unsigned b = getBorder(e);
			if (b > drb)
				drb = b;
			TPair pair(record_num, b);
			appendValue(dep, pair);
		}

		else // end of dependent region
		{
			//deleteEntries(records, recs_to_delete);
			return (i);
		}
	}

	//deleteEntries(records, recs_to_delete);
	return length(records);
}

// ----------------------------------------------------------------------------
// Function processDeltaEventsOnReference()
// ----------------------------------------------------------------------------
//main
template<typename TDeltaEvents, typename TSequence>
int processDeltaEventsOnReference(TDeltaEvents & records, TSequence & ref)
{
	String<unsigned> recs_to_delete;
	sort(records, CompareByPosAndTypeLessThan_()); // sort by reference position(first) and type(second) and value(third)

	unsigned start = 0;
	unsigned end = length(records);
	while(start < end)
	{
		String<unsigned> tmp_recs_to_delete;



		DependentRegion dr;
		unsigned new_start = getNextDR(dr, tmp_recs_to_delete, records, start);
		int offset = processDR(dr, ref); // inside here, no record is deleted

		// update records from start to new_start todo:: why does it not work through references???
		unsigned skip = 0;
		for(unsigned i = start; i < new_start; ++i)
		{
			if (!(isIn(tmp_recs_to_delete, i)))
				assignValue(records, i, dr.records[i-start-skip]);
			else
				++skip;
		}
		append(recs_to_delete, tmp_recs_to_delete);

		if(offset!=0)
			applyOffset(records, offset, new_start);
		start = new_start;
		end = length(records); // update because number of records might change!
	}
	deleteEntries(records, recs_to_delete);
	return 0;
}
#endif /* APPS_LAGAN_PROCESSEVENTS_H_ */
