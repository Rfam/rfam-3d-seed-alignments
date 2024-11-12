import React, { useState, useMemo, useEffect } from "react";
import { ChevronDown, ChevronRight, Search } from "lucide-react";
import { fetchRfamData } from "../services/dataService";

const RfamDashboard = () => {
  const [searchTerm, setSearchTerm] = useState("");
  const [data, setData] = useState([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);
  const [expandedStructures, setExpandedStructures] = useState({});
  const [filters, setFilters] = useState({
    hasBasePairs: false,
    hasPseudoknots: false,
    chainStatus: "all",
    maxResolution: 4.0,
  });

  useEffect(() => {
    const loadData = async () => {
      try {
        setLoading(true);
        const rfamData = await fetchRfamData();
        setData(rfamData);
      } catch (err) {
        setError(err.message);
      } finally {
        setLoading(false);
      }
    };

    loadData();
  }, []);

  const toggleStructure = (familyId, pdbId) => {
    setExpandedStructures((prev) => ({
      ...prev,
      [`${familyId}-${pdbId}`]: !prev[`${familyId}-${pdbId}`],
    }));
  };

  const summary = useMemo(() => {
    const stats = {
      totalFamilies: data.length,
      familiesWithNewStructures: 0,
      familiesWithNewSequences: 0,
      familiesWithHighResolution: 0,
      statusBreakdown: {
        complete: 0,
        curatable: 0,
        skipped: 0,
        no_matches: 0,
      },
    };

    data.forEach((family) => {
      let hasNewStructure = false;
      let hasNewSequence = false;
      let hasHighResolution = false;

      // Increment status count
      stats.statusBreakdown[family.status] =
        (stats.statusBreakdown[family.status] || 0) + 1;

      family.structures.forEach((structure) => {
        if (structure.resolution <= 4.0) {
          hasHighResolution = true;
        }
        structure.chains.forEach((chain) => {
          if (chain.status === "new_structure") hasNewStructure = true;
          if (chain.status === "new_sequence") hasNewSequence = true;
        });
      });

      if (hasNewStructure) stats.familiesWithNewStructures++;
      if (hasNewSequence) stats.familiesWithNewSequences++;
      if (hasHighResolution) stats.familiesWithHighResolution++;
    });

    return stats;
  }, [data]);

  const getStatusColor = (status) => {
    switch (status) {
      // Family statuses
      case "complete":
        return "bg-green-100 text-green-800";
      case "curatable":
        return "bg-blue-100 text-blue-800";
      case "skipped":
        return "bg-yellow-100 text-yellow-800";
      case "no_matches":
        return "bg-gray-100 text-gray-800";
      // Chain statuses
      case "already_present":
        return "bg-green-100 text-green-800";
      case "new_sequence":
        return "bg-blue-100 text-blue-800";
      case "new_structure":
        return "bg-purple-100 text-purple-800";
      default:
        return "bg-gray-100 text-gray-800";
    }
  };

  const getStatusLabel = (status) => {
    switch (status) {
      case "no_matches":
        return "No Matches";
      default:
        return status
          .split("_")
          .map((word) => word.charAt(0).toUpperCase() + word.slice(1))
          .join(" ");
    }
  };

  const filteredData = useMemo(() => {
    return data.filter((family) => {
      const matchesFamily =
        family.familyId.toLowerCase().includes(searchTerm.toLowerCase()) ||
        family.familyName.toLowerCase().includes(searchTerm.toLowerCase());

      const matchesPDB = family.structures.some((structure) =>
        structure.pdbId.toLowerCase().includes(searchTerm.toLowerCase()),
      );

      const matchesSearch = matchesFamily || matchesPDB;

      const matchesFilters = family.structures.some((structure) => {
        const meetsResolution = structure.resolution <= filters.maxResolution;

        const matchingChains = structure.chains.some((chain) => {
          const matchesBasePairs = !filters.hasBasePairs || chain.hasBasePairs;
          const matchesPseudoknots =
            !filters.hasPseudoknots || chain.hasPseudoknots;
          const matchesStatus =
            filters.chainStatus === "all" ||
            chain.status === filters.chainStatus;

          return matchesBasePairs && matchesPseudoknots && matchesStatus;
        });

        return meetsResolution && matchingChains;
      });

      return (
        matchesSearch && (family.structures.length === 0 || matchesFilters)
      );
    });
  }, [data, searchTerm, filters]);

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
              <div className={`w-3 h-3 rounded-full bg-green-100`}></div>
              <span>Complete: {summary.statusBreakdown.complete}</span>
            </div>
            <div className="flex items-center gap-2">
              <div className={`w-3 h-3 rounded-full bg-blue-100`}></div>
              <span>Curatable: {summary.statusBreakdown.curatable}</span>
            </div>
            <div className="flex items-center gap-2">
              <div className={`w-3 h-3 rounded-full bg-yellow-100`}></div>
              <span>Skipped: {summary.statusBreakdown.skipped}</span>
            </div>
            <div className="flex items-center gap-2">
              <div className={`w-3 h-3 rounded-full bg-gray-100`}></div>
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

          <div className="flex flex-wrap gap-4">
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

            <select
              value={filters.chainStatus}
              onChange={(e) =>
                setFilters((prev) => ({ ...prev, chainStatus: e.target.value }))
              }
              className="border rounded-md px-2 py-1"
            >
              <option value="all">All Chain Statuses</option>
              <option value="already_present">Already Present</option>
              <option value="new_sequence">New Sequence</option>
              <option value="new_structure">New Structure</option>
            </select>

            <label className="flex items-center gap-2">
              <span>Max Resolution (Å):</span>
              <input
                type="number"
                value={filters.maxResolution}
                onChange={(e) =>
                  setFilters((prev) => ({
                    ...prev,
                    maxResolution: parseFloat(e.target.value) || 4.0,
                  }))
                }
                className="border rounded-md px-2 py-1 w-20"
                step="0.1"
              />
            </label>
          </div>
        </div>
      </div>

      <div className="space-y-4">
        {filteredData.map((family) => (
          <div
            key={family.familyId}
            className="border rounded-lg p-4 bg-white shadow-sm"
          >
            <div className="flex justify-between items-center mb-2">
              <h3 className="text-lg font-medium">
                {family.familyId} - {family.familyName}
              </h3>
              <span
                className={`px-3 py-1 rounded-full text-sm ${getStatusColor(family.status)}`}
              >
                {getStatusLabel(family.status)}
              </span>
            </div>

            {family.structures.length > 0 && (
              <div className="mt-2">
                <h4 className="text-sm font-medium text-gray-500 mb-2">
                  Structures
                </h4>
                <div className="space-y-2">
                  {family.structures.map((structure) => (
                    <div key={structure.pdbId} className="border rounded p-2">
                      <div
                        className="flex items-center cursor-pointer"
                        onClick={() =>
                          toggleStructure(family.familyId, structure.pdbId)
                        }
                      >
                        {expandedStructures[
                          `${family.familyId}-${structure.pdbId}`
                        ] ? (
                          <ChevronDown className="h-4 w-4 mr-2" />
                        ) : (
                          <ChevronRight className="h-4 w-4 mr-2" />
                        )}
                        <span className="font-medium flex-grow">
                          PDB: {structure.pdbId}
                        </span>
                        <a
                          href={`https://www.rcsb.org/structure/${structure.pdbId}`}
                          target="_blank"
                          rel="noopener noreferrer"
                          className="text-blue-600 hover:text-blue-800 mr-4"
                          onClick={(e) => e.stopPropagation()}
                        >
                          View in RCSB
                        </a>
                        <span className="text-sm text-gray-500 mr-4">
                          {structure.method}
                        </span>
                        <span className="text-sm text-gray-500 mr-4">
                          Resolution: {structure.resolution}Å
                        </span>
                        <span className="text-sm text-gray-500">
                          {structure.chains.length} chain
                          {structure.chains.length !== 1 ? "s" : ""}
                        </span>
                      </div>

                      {expandedStructures[
                        `${family.familyId}-${structure.pdbId}`
                      ] && (
                        <div className="mt-2 ml-6">
                          <div className="space-y-2">
                            {structure.chains.map((chain) => (
                              <div
                                key={chain.chainId}
                                className="border rounded p-3"
                              >
                                <div className="flex flex-wrap gap-2 items-center">
                                  <span className="font-medium min-w-20">
                                    Chain {chain.chainId}
                                  </span>
                                  <span
                                    className={`px-2 py-0.5 rounded-full text-xs ${getStatusColor(chain.status)}`}
                                  >
                                    {getStatusLabel(chain.status)}
                                  </span>
                                  <div className="flex gap-4 ml-auto">
                                    <span
                                      className={`text-sm ${chain.hasBasePairs ? "text-green-600" : "text-gray-400"}`}
                                    >
                                      {chain.hasBasePairs ? "✓" : "✗"} Base
                                      pairs
                                    </span>
                                    <span
                                      className={`text-sm ${chain.hasPseudoknots ? "text-green-600" : "text-gray-400"}`}
                                    >
                                      {chain.hasPseudoknots ? "✓" : "✗"}{" "}
                                      Pseudoknots
                                    </span>
                                  </div>
                                </div>
                              </div>
                            ))}
                          </div>
                        </div>
                      )}
                    </div>
                  ))}
                </div>
              </div>
            )}
          </div>
        ))}
      </div>
    </div>
  );
};

export default RfamDashboard;
