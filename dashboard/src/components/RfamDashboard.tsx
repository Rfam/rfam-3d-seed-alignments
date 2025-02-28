// RfamDashboard.tsx - Updated without summary information
import React, { useState, useMemo, useEffect } from "react";
import { Search } from "lucide-react";
import { RfamService } from "../services/dataService";
import { RfamFamily, RfamFilters, RfamSequence } from "../types";
import CompactFamilyCard from "./CompactFamilyCard";

const fetchRfamData = async (): Promise<RfamFamily[]> => {
  console.log("[Dashboard] Starting to fetch Rfam data");
  const rfamService = RfamService.getInstance();
  try {
    const data = await rfamService.getData();
    console.log("[Dashboard] Successfully fetched data, items:", data.length);
    return data;
  } catch (error) {
    console.error("[Dashboard] Error fetching Rfam data:", error);
    throw error;
  }
};

const RfamDashboard: React.FC = () => {
  const [searchTerm, setSearchTerm] = useState<string>("");
  const [data, setData] = useState<RfamFamily[]>([]);
  const [loading, setLoading] = useState<boolean>(true);
  const [error, setError] = useState<string | null>(null);
  // No filters except search

  useEffect(() => {
    const loadData = async () => {
      console.log("[Dashboard] Starting to load data");
      setLoading(true);
      try {
        console.log("[Dashboard] Calling fetchRfamData()");
        const result = await fetchRfamData();
        console.log("[Dashboard] Setting data state with result");
        setData(result);
      } catch (err) {
        console.error("[Dashboard] Error in loadData:", err);
        setError(err instanceof Error ? err.message : "Failed to load data");
      } finally {
        console.log("[Dashboard] Finished loading data, setting loading=false");
        setLoading(false);
      }
    };
    loadData();
  }, []);

  // Helper function to check if a family matches the search term
  const familyMatchesSearch = (family: RfamFamily): boolean => {
    // Return all data if no search term
    if (!searchTerm) return true;

    const lowercasedSearch = searchTerm.toLowerCase();

    // Check if family ID or name matches
    if (
      family.familyId.toLowerCase().includes(lowercasedSearch) ||
      family.familyName.toLowerCase().includes(lowercasedSearch)
    ) {
      return true;
    }

    // Check if any structure PDB ID matches
    return family.sequences.some((sequence) =>
      sequence.structures.some((structure) =>
        structure.pdbId.toLowerCase().includes(lowercasedSearch),
      ),
    );
  };

  const filteredData = useMemo<RfamFamily[]>(() => {
    return data.filter(familyMatchesSearch);
  }, [data, searchTerm]);

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
      <div className="mb-6">
        <h1 className="text-2xl font-bold mb-6">Rfam Analysis Dashboard</h1>

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

        {/* All filters removed except search */}

        {/* Results Count */}
        <div className="mt-6 mb-2 text-sm text-gray-600">
          Showing {filteredData.length}{" "}
          {filteredData.length === 1 ? "family" : "families"}
        </div>

        {/* Family List */}
        <div className="space-y-2">
          {filteredData.length > 0 ? (
            filteredData.map((family) => (
              <CompactFamilyCard key={family.familyId} family={family} />
            ))
          ) : (
            <div className="bg-gray-50 p-8 text-center rounded-lg">
              <p className="text-gray-500">
                No families match the current filters
              </p>
            </div>
          )}
        </div>
      </div>
    </div>
  );
};

export default RfamDashboard;
