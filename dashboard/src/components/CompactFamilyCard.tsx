import React, { useState } from "react";
import { Info, ChevronDown, ChevronRight, ExternalLink } from "lucide-react";
import { RfamFamily, RfamSequence, Structure, Chain } from "../types";

interface CompactFamilyCardProps {
  family: RfamFamily;
}

const CompactFamilyCard: React.FC<CompactFamilyCardProps> = ({ family }) => {
  const [isExpanded, setIsExpanded] = useState(false);

  // Calculate summary statistics from sequences with safe guards for missing data
  const stats = {
    totalStructures: family.metadata?.totalStructures || 0,
    newStructures: (family.sequences || []).filter(
      (seq) => seq && seq.status === "new_structure",
    ).length,
    // Count both new_structure and new_sequence for newSequences
    newSequences: (family.sequences || []).filter(
      (seq) =>
        seq &&
        (seq.status === "new_structure" || seq.status === "new_sequence"),
    ).length,
    structuresWithBasePairs: (family.sequences || []).filter(
      (seq) => seq && seq.metadata && seq.metadata.hasBasePairs,
    ).length,
    structuresWithPseudoknots: (family.sequences || []).filter(
      (seq) => seq && seq.metadata && seq.metadata.hasPseudoknots,
    ).length,
  };

  const getStatusColor = (status: string): string => {
    const colors: Record<string, string> = {
      complete: "bg-green-100 text-green-800",
      curatable: "bg-blue-100 text-blue-800",
      skipped: "bg-yellow-100 text-yellow-800",
      no_matches: "bg-gray-100 text-gray-800",
      new_structure: "bg-purple-100 text-purple-800",
      new_sequence: "bg-blue-100 text-blue-800",
      committed: "bg-gray-100 text-gray-800",
      error: "bg-red-100 text-red-800",
      unknown: "bg-gray-100 text-gray-800",
    };
    return colors[status] || "bg-gray-100 text-gray-800";
  };

  const shouldShowChainStatus = (status: string | undefined): boolean => {
    if (!status) return false;
    return !["new_structure", "new_sequence"].includes(status);
  };

  const hasExpandableContent = stats.totalStructures > 0;

  return (
    <div className="border rounded-lg p-3 bg-white shadow-sm">
      <div
        className={`flex flex-wrap gap-2 items-start ${hasExpandableContent ? "cursor-pointer" : ""}`}
        onClick={() => hasExpandableContent && setIsExpanded(!isExpanded)}
      >
        {/* Left column - Basic info */}
        <div className="flex-grow min-w-[200px] max-w-[300px]">
          <div className="flex items-center gap-2">
            {hasExpandableContent &&
              (isExpanded ? (
                <ChevronDown className="h-4 w-4 text-gray-500" />
              ) : (
                <ChevronRight className="h-4 w-4 text-gray-500" />
              ))}
            <h3
              className={`font-medium text-gray-900 ${!hasExpandableContent ? "ml-6" : ""}`}
            >
              {family.familyId}
            </h3>
            <span
              className={`px-2 py-0.5 text-xs rounded-full ${getStatusColor(family.status)}`}
            >
              {family.status}
            </span>
          </div>
          <p
            className="text-sm text-gray-600 truncate ml-6"
            title={family.familyName}
          >
            {family.familyName}
          </p>
        </div>

        {/* Middle column removed */}

        {/* Right column - Feature indicators */}
        {stats.totalStructures > 0 && (
          <div className="flex gap-3 items-center ml-auto">
            <div className="flex flex-col items-center gap-1">
              <div
                className={`h-2 w-2 rounded-full ${stats.structuresWithBasePairs > 0 ? "bg-green-500" : "bg-gray-200"}`}
              />
              <span className="text-xs text-gray-600">Base Pairs</span>
            </div>
            <div className="flex flex-col items-center gap-1">
              <div
                className={`h-2 w-2 rounded-full ${stats.structuresWithPseudoknots > 0 ? "bg-green-500" : "bg-gray-200"}`}
              />
              <span className="text-xs text-gray-600">Pseudoknots</span>
            </div>
            <a
              href={`https://rfam.org/family/${family.familyId}`}
              target="_blank"
              rel="noopener noreferrer"
              className="ml-2 p-1 text-blue-600 hover:text-blue-800"
              title="View family in Rfam database"
              onClick={(e) => e.stopPropagation()}
            >
              <Info className="h-4 w-4" />
            </a>
            {family.status !== "skipped" && (
              <a
                href={`https://ftp.ebi.ac.uk/pub/databases/Rfam/.preview/3d/alignments/${family.familyId}.updated.sto`}
                target="_blank"
                rel="noopener noreferrer"
                className="p-1 text-blue-600 hover:text-blue-800"
                title="View family alignment"
                onClick={(e) => e.stopPropagation()}
              >
                <svg
                  className="h-4 w-4"
                  xmlns="http://www.w3.org/2000/svg"
                  viewBox="0 0 24 24"
                  fill="none"
                  stroke="currentColor"
                  strokeWidth="2"
                  strokeLinecap="round"
                  strokeLinejoin="round"
                >
                  <line x1="21" y1="10" x2="3" y2="10"></line>
                  <line x1="21" y1="6" x2="3" y2="6"></line>
                  <line x1="21" y1="14" x2="3" y2="14"></line>
                  <line x1="21" y1="18" x2="3" y2="18"></line>
                </svg>
              </a>
            )}
          </div>
        )}
      </div>

      {/* Expanded View - Sequences */}
      {isExpanded && hasExpandableContent && (
        <div className="mt-4 ml-6 space-y-4">
          {family.sequences && family.sequences.length > 0 ? (
            family.sequences.map((sequence, seqIndex) =>
              sequence ? (
                <div
                  key={`${sequence.sequenceId || seqIndex}`}
                  className="border rounded-lg p-3"
                >
                  <div className="flex items-center justify-between mb-2">
                    <div className="font-medium text-sm text-gray-700">
                      {sequence.sequenceId || "Unknown Sequence"}
                    </div>
                    <span
                      className={`px-2 py-0.5 text-xs rounded-full ${getStatusColor(sequence.status || "unknown")}`}
                    >
                      {(sequence.status || "unknown")
                        .split("_")
                        .map(
                          (word) =>
                            word.charAt(0).toUpperCase() + word.slice(1),
                        )
                        .join(" ")}
                    </span>
                  </div>
                  <div className="space-y-2">
                    {sequence.structures && sequence.structures.length > 0 ? (
                      sequence.structures.map((structure) =>
                        structure.chains && structure.chains.length > 0 ? (
                          structure.chains.map((chain, chainIndex) => (
                            <div
                              key={`${structure.pdbId || "unknown"}-${chain.chainId || "unknown"}-${chainIndex}`}
                              className="flex items-center gap-4 text-sm pl-2"
                            >
                              {chain.status &&
                                shouldShowChainStatus(chain.status) && (
                                  <span
                                    className={`px-2 py-0.5 rounded-full text-xs ${getStatusColor(chain.status)}`}
                                  >
                                    {chain.status
                                      .split("_")
                                      .map(
                                        (word) =>
                                          word.charAt(0).toUpperCase() +
                                          word.slice(1),
                                      )
                                      .join(" ")}
                                  </span>
                                )}
                              <a
                                href={`https://www.rcsb.org/structure/${structure.pdbId || "unknown"}`}
                                target="_blank"
                                rel="noopener noreferrer"
                                className="flex items-center gap-1 text-blue-600 hover:text-blue-800"
                                onClick={(e) => e.stopPropagation()}
                              >
                                {structure.pdbId || "unknown"}_
                                {chain.chainId || "unknown"}
                                <ExternalLink className="h-3 w-3" />
                              </a>
                              <span className="text-gray-500">
                                {structure.resolution &&
                                structure.resolution > 0
                                  ? `${structure.resolution.toFixed(1)}Å`
                                  : "N/A"}
                              </span>
                              <div className="flex gap-3 ml-auto">
                                <span
                                  className={
                                    chain.hasBasePairs
                                      ? "text-green-600"
                                      : "text-gray-400"
                                  }
                                >
                                  {chain.hasBasePairs ? "✓" : "✗"} Base Pairs
                                </span>
                                <span
                                  className={
                                    chain.hasPseudoknots
                                      ? "text-green-600"
                                      : "text-gray-400"
                                  }
                                >
                                  {chain.hasPseudoknots ? "✓" : "✗"} Pseudoknots
                                </span>
                              </div>
                            </div>
                          ))
                        ) : (
                          <div
                            key={`no-chains-${structure.pdbId || "unknown"}`}
                            className="text-gray-500 text-sm italic pl-2"
                          >
                            No chains available for structure{" "}
                            {structure.pdbId || "unknown"}
                          </div>
                        ),
                      )
                    ) : (
                      <div className="text-gray-500 text-sm italic pl-2">
                        No structures available
                      </div>
                    )}
                  </div>
                </div>
              ) : (
                <div
                  key={`missing-sequence-${seqIndex}`}
                  className="border rounded-lg p-3"
                >
                  <div className="text-gray-500 italic">
                    Missing sequence data
                  </div>
                </div>
              ),
            )
          ) : (
            <div className="text-gray-500 italic p-3">
              No sequence data available
            </div>
          )}
        </div>
      )}
    </div>
  );
};

export default CompactFamilyCard;
