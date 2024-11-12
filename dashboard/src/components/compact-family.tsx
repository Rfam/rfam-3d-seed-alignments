import React, { useState } from "react";
import { Info, ChevronDown, ChevronRight, ExternalLink } from "lucide-react";
import type { RfamFamily, Structure, Chain, RfamSequence } from "../types";

interface CompactFamilyCardProps {
  family: RfamFamily;
}

const CompactFamilyCard: React.FC<CompactFamilyCardProps> = ({ family }) => {
  const [isExpanded, setIsExpanded] = useState(false);

  // Calculate summary statistics from sequences
  const stats = {
    totalStructures: family.metadata.totalStructures,
    newStructures: family.sequences.filter(
      (seq) => seq.status === "new_structure",
    ).length,
    newSequences: family.sequences.filter(
      (seq) => seq.status === "new_structure" || seq.status === "new_sequence",
    ).length,
    structuresWithBasePairs: family.sequences.filter(
      (seq) => seq.metadata.hasBasePairs,
    ).length,
    structuresWithPseudoknots: family.sequences.filter(
      (seq) => seq.metadata.hasPseudoknots,
    ).length,
  };

  const getStatusColor = (
    status: RfamFamily["status"] | Chain["status"] | RfamSequence["status"],
  ): string => {
    const colors: Record<string, string> = {
      complete: "bg-green-100 text-green-800",
      incomplete: "bg-yellow-100 text-yellow-800",
      curate: "bg-blue-100 text-blue-800",
      skipped: "bg-gray-100 text-gray-800",
      no_matches: "bg-gray-100 text-gray-800",
      new_structure: "bg-purple-100 text-purple-800",
      new_sequence: "bg-blue-100 text-blue-800",
      already_present: "bg-gray-100 text-gray-800",
      error: "bg-red-100 text-red-800",
    };
    return colors[status] || "bg-gray-100 text-gray-800";
  };

  const shouldShowChainStatus = (status: Chain["status"]): boolean => {
    return !["new_structure", "new_sequence"].includes(status);
  };

  const formatResolution = (resolution: number | null): string => {
    if (resolution === null) return "N/A";
    return `${resolution.toFixed(1)}Å`;
  };

  const formatMethod = (method: string | null): string => {
    return method || "N/A";
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
            <a
              href={`https://rfam.org/family/${family.familyId}`}
              target="_blank"
              rel="noopener noreferrer"
              className="text-blue-600 hover:text-blue-800"
              title="View family in Rfam"
              onClick={(e) => e.stopPropagation()}
            >
              <Info className="h-4 w-4" />
            </a>
          </div>
          <p
            className="text-sm text-gray-600 truncate ml-6"
            title={family.familyName}
          >
            {family.familyName}
          </p>
        </div>

        {/* Middle column - Structure summary */}
        <div className="flex-grow min-w-[200px]">
          {stats.totalStructures > 0 ? (
            <div className="grid grid-cols-2 gap-x-4 gap-y-1 text-sm">
              <div className="flex items-center gap-1">
                <span className="text-gray-600">Structures:</span>
                <span className="font-medium">{stats.totalStructures}</span>
              </div>
              <div className="flex items-center gap-1">
                <span className="text-gray-600">New Struct:</span>
                <span
                  className={`font-medium ${stats.newStructures > 0 ? "text-purple-600" : ""}`}
                >
                  {stats.newStructures}
                </span>
              </div>
              <div className="flex items-center gap-1">
                <span className="text-gray-600">New Seq:</span>
                <span
                  className={`font-medium ${stats.newSequences > 0 ? "text-blue-600" : ""}`}
                >
                  {stats.newSequences}
                </span>
              </div>
            </div>
          ) : (
            <span className="text-sm text-gray-500">
              No structures available
            </span>
          )}
        </div>

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
          </div>
        )}
      </div>

      {/* Expanded View - Sequences */}
      {isExpanded && hasExpandableContent && (
        <div className="mt-4 ml-6 space-y-4">
          {family.sequences.map((sequence) => (
            <div key={sequence.sequenceId} className="border rounded-lg p-3">
              <div className="flex items-center justify-between mb-2">
                <div className="font-medium text-sm text-gray-700">
                  {sequence.sequenceId}
                  <span className="text-gray-500 ml-2">
                    ({sequence.structures.length} structure
                    {sequence.structures.length !== 1 ? "s" : ""})
                  </span>
                </div>
                <span
                  className={`px-2 py-0.5 text-xs rounded-full ${getStatusColor(sequence.status)}`}
                >
                  {sequence.status
                    .split("_")
                    .map((word) => word.charAt(0).toUpperCase() + word.slice(1))
                    .join(" ")}
                </span>
              </div>
              <div className="space-y-2">
                {sequence.structures.map((structure) =>
                  structure.chains.map((chain) => (
                    <div
                      key={`${structure.pdbId}-${chain.chainId}`}
                      className="flex items-center gap-4 text-sm pl-2"
                    >
                      {shouldShowChainStatus(chain.status) && (
                        <span
                          className={`px-2 py-0.5 rounded-full text-xs ${getStatusColor(chain.status)}`}
                        >
                          {chain.status
                            .split("_")
                            .map(
                              (word) =>
                                word.charAt(0).toUpperCase() + word.slice(1),
                            )
                            .join(" ")}
                        </span>
                      )}
                      <a
                        href={`https://www.rcsb.org/structure/${structure.pdbId}`}
                        target="_blank"
                        rel="noopener noreferrer"
                        className="flex items-center gap-1 text-blue-600 hover:text-blue-800"
                        onClick={(e) => e.stopPropagation()}
                      >
                        {structure.pdbId}:{chain.chainId}
                        <ExternalLink className="h-3 w-3" />
                      </a>
                      <div className="flex items-center gap-4">
                        <span className="text-gray-500">
                          {formatResolution(structure.resolution)}
                        </span>
                        {structure.method && (
                          <span className="text-gray-500">
                            {formatMethod(structure.method)}
                          </span>
                        )}
                      </div>
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
                  )),
                )}
              </div>
            </div>
          ))}
        </div>
      )}
    </div>
  );
};

export default CompactFamilyCard;
