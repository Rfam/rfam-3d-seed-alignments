import React, { useState, useMemo, useEffect } from "react";
import { Search } from "lucide-react";
import { RfamService } from "../services/dataService";
import type {
  RfamFamily,
  RfamFilters,
  RfamSummary,
  RfamSequence,
} from "../types";
import CompactFamilyCard from "./compactFamilyCard";

const fetchRfamData = async (): Promise<RfamFamily[]> => {
  const rfamService = RfamService.getInstance();
  return rfamService.getData();
};

const RfamDashboard: React.FC = () => {
  const [searchTerm, setSearchTerm] = useState<string>("");
  const [data, setData] = useState<RfamFamily[]>([]);
  const [loading, setLoading] = useState<boolean>(true);
  const [error, setError] = useState<string | null>(null);
  const [filters, setFilters] = useState<RfamFilters>({
    hasBasePairs: false,
    hasPseudoknots: false,
    familyStatus: "all",
    sequenceStatus: "all",
    chainStatus: "all",
  });

  useEffect(() => {
    const loadData = async () => {
      try {
        const result = await fetchRfamData();
        setData(result);
      } catch (err) {
        setError(err instanceof Error ? err.message : "Failed to load data");
      } finally {
        setLoading(false);
      }
    };
    loadData();
  }, []);

  const hasActiveStructuralFilters = (): boolean => {
    return filters.hasBasePairs || filters.hasPseudoknots;
  };

  // Helper function to check if a sequence matches the filters
  const sequenceMatchesFilters = (sequence: RfamSequence): boolean => {
    const matchesStatus =
      filters.sequenceStatus === "all" ||
      sequence.status === filters.sequenceStatus;
    const matchesBasePairs =
      !filters.hasBasePairs || sequence.metadata.hasBasePairs;
    const matchesPseudoknots =
      !filters.hasPseudoknots || sequence.metadata.hasPseudoknots;

    return matchesStatus && matchesBasePairs && matchesPseudoknots;
  };

  // Helper function to check if a family matches current filters
  const familyMatchesFilters = (family: RfamFamily): boolean => {
    // Family status filter
    const matchesFamilyStatus =
      filters.familyStatus === "all" || family.status === filters.familyStatus;

    // Search term filter
    const matchesSearch =
      !searchTerm ||
      family.familyId.toLowerCase().includes(searchTerm.toLowerCase()) ||
      family.familyName.toLowerCase().includes(searchTerm.toLowerCase()) ||
      family.sequences.some((sequence) =>
        sequence.structures.some((structure) =>
          structure.pdbId.toLowerCase().includes(searchTerm.toLowerCase()),
        ),
      );

    // For structural filters, we need sequence matching
    if (hasActiveStructuralFilters()) {
      const hasMatchingSequence = family.sequences.some(sequenceMatchesFilters);
      return matchesFamilyStatus && matchesSearch && hasMatchingSequence;
    }

    // If no structural filters are active, just check family status and search
    return matchesFamilyStatus && matchesSearch;
  };

  const filteredData = useMemo<RfamFamily[]>(() => {
    return data.filter(familyMatchesFilters);
  }, [data, searchTerm, filters]);

  const summary = useMemo<RfamSummary>(() => {
    const stats: RfamSummary = {
      totalFamilies: filteredData.length,
      totalSequences: 0,
      totalStructures: 0,
      familiesWithNewStructures: 0,
      familiesWithNewSequences: 0,
      familiesWithHighResolution: 0,
      statusBreakdown: {
        complete: 0,
        curatable: 0,
        skipped: 0,
        no_matches: 0,
      },
      sequenceBreakdown: {
        already_present: 0,
        new_sequence: 0,
        new_structure: 0,
      },
    };

    filteredData.forEach((family) => {
      stats.statusBreakdown[family.status]++;
      stats.totalSequences += family.sequences.length;
      stats.totalStructures += family.metadata.totalStructures;

      const hasNewStructure = family.sequences.some(
        (seq) => seq.status === "new_structure",
      );
      const hasNewSequence = family.sequences.some(
        (seq) => seq.status === "new_sequence",
      );
      const hasHighResolution = family.metadata.bestResolution <= 4.0;

      if (hasNewStructure) stats.familiesWithNewStructures++;
      if (hasNewSequence) stats.familiesWithNewSequences++;
      if (hasHighResolution) stats.familiesWithHighResolution++;

      family.sequences.forEach((sequence) => {
        stats.sequenceBreakdown[sequence.status]++;
      });
    });

    return stats;
  }, [filteredData]);

  if (loading) {
    return (
      <div className="max-w-4xl mx-auto p-6">
        <div className="flex items-center justify-center h-64">
          <div className="text-gray-600">Loading data...</div>
        </div>
      </div>
    );
  }

  if (error) {
    return (
      <div className="max-w-4xl mx-auto p-6">
        <div className="bg-red-50 border border-red-200 rounded-lg p-4">
          <h2 className="text-red-800 font-medium">Error Loading Data</h2>
          <p className="text-red-600">{error}</p>
          <button
            onClick={() => window.location.reload()}
            className="mt-2 text-red-700 hover:text-red-800 underline"
          >
            Retry
          </button>
        </div>
      </div>
    );
  }

  return (
    <div className="max-w-4xl mx-auto p-6">
      <div className="mb-8">
        <h1 className="text-2xl font-bold mb-6">
          Rfam Analysis Status Dashboard
        </h1>

        {/* Main Summary Stats */}
        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-4 mb-6">
          <div className="bg-blue-50 p-4 rounded-lg">
            <div className="text-2xl font-bold">
              {summary.familiesWithNewSequences}
            </div>
            <div className="text-sm text-gray-600">
              Families with New Sequences
            </div>
          </div>
          <div className="bg-purple-50 p-4 rounded-lg">
            <div className="text-2xl font-bold">
              {summary.familiesWithNewStructures}
            </div>
            <div className="text-sm text-gray-600">
              Families with New Structures
            </div>
          </div>
          <div className="bg-green-50 p-4 rounded-lg">
            <div className="text-2xl font-bold">
              {summary.familiesWithHighResolution}
            </div>
            <div className="text-sm text-gray-600">
              Families with High Resolution (≤4.0Å)
            </div>
          </div>
          <div className="bg-gray-50 p-4 rounded-lg">
            <div className="text-2xl font-bold">{summary.totalFamilies}</div>
            <div className="text-sm text-gray-600">Total Families</div>
          </div>
        </div>

        {/* Status Breakdown */}
        <div className="bg-white rounded-lg shadow-sm p-4 mb-6">
          <h2 className="text-lg font-medium mb-4">Family Status Breakdown</h2>
          <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
            <div className="flex items-center gap-2">
              <div className="w-3 h-3 rounded-full bg-green-100"></div>
              <span>Complete: {summary.statusBreakdown.complete}</span>
            </div>
            <div className="flex items-center gap-2">
              <div className="w-3 h-3 rounded-full bg-blue-100"></div>
              <span>Curatable: {summary.statusBreakdown.curatable}</span>
            </div>
            <div className="flex items-center gap-2">
              <div className="w-3 h-3 rounded-full bg-yellow-100"></div>
              <span>Skipped: {summary.statusBreakdown.skipped}</span>
            </div>
            <div className="flex items-center gap-2">
              <div className="w-3 h-3 rounded-full bg-gray-100"></div>
              <span>No Matches: {summary.statusBreakdown.no_matches}</span>
            </div>
          </div>
        </div>

        {/* Search and Filters Section */}
        <div className="space-y-4">
          <div className="relative">
            <Search className="absolute left-2 top-2.5 h-4 w-4 text-gray-500" />
            <input
              type="text"
              placeholder="Search by family ID, name, or PDB ID..."
              className="pl-8 pr-4 py-2 w-full border rounded-md"
              value={searchTerm}
              onChange={(e) => setSearchTerm(e.target.value)}
            />
          </div>
        </div>

        {/* Filters section with active filter indicator */}
        <div className="space-y-4">
          <div className="flex flex-wrap gap-4">
            <div className="w-full">
              {hasActiveStructuralFilters() && (
                <div className="mb-4 p-2 bg-yellow-50 border border-yellow-200 rounded-md text-yellow-800">
                  Note: Complete and Skipped families are hidden when structural
                  filters are active
                </div>
              )}
            </div>

            {/* Base Pairs Filter */}
            <label className="flex items-center gap-2">
              <input
                type="checkbox"
                checked={filters.hasBasePairs}
                onChange={(e) =>
                  setFilters((prev) => ({
                    ...prev,
                    hasBasePairs: e.target.checked,
                  }))
                }
                className="rounded border-gray-300"
              />
              <span>Has Base Pairs</span>
            </label>

            {/* Pseudoknots Filter */}
            <label className="flex items-center gap-2">
              <input
                type="checkbox"
                checked={filters.hasPseudoknots}
                onChange={(e) =>
                  setFilters((prev) => ({
                    ...prev,
                    hasPseudoknots: e.target.checked,
                  }))
                }
                className="rounded border-gray-300"
              />
              <span>Has Pseudoknots</span>
            </label>
          </div>

          {/* Status Filters */}
          <div className="flex gap-4 w-full">
            <select
              value={filters.familyStatus}
              onChange={(e) =>
                setFilters((prev) => ({
                  ...prev,
                  familyStatus: e.target.value as RfamFamily["status"] | "all",
                }))
              }
              className="border rounded-md px-2 py-1"
            >
              <option value="all">All Family Statuses</option>
              <option value="complete">Complete</option>
              <option value="incomplete">Incomplete</option>
              <option value="curate">Curate</option>
              <option value="skipped">Skipped</option>
              <option value="no_matches">No Matches</option>
            </select>

            <select
              value={filters.sequenceStatus}
              onChange={(e) =>
                setFilters((prev) => ({
                  ...prev,
                  sequenceStatus: e.target.value as
                    | RfamSequence["status"]
                    | "all",
                }))
              }
              className="border rounded-md px-2 py-1"
            >
              <option value="all">All Sequence Statuses</option>
              <option value="already_present">Already Present</option>
              <option value="new_sequence">New Sequence</option>
              <option value="new_structure">New Structure</option>
            </select>
          </div>
        </div>

        {/* Family List */}
        <div className="space-y-2">
          {filteredData.map((family) => (
            <CompactFamilyCard key={family.familyId} family={family} />
          ))}
        </div>
      </div>
    </div>
  );
};

export default RfamDashboard;
