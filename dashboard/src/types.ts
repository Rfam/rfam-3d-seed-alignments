// Chain interface for individual RNA chain data
export interface Chain {
  chainId: string;
  status:
    | "already_present"
    | "new_sequence"
    | "new_structure"
    | "error"
    | "skipped";
  hasBasePairs: boolean;
  hasPseudoknots: boolean;
}

// Structure interface for PDB structure data
export interface Structure {
  pdbId: string;
  method: string | null;
  resolution: number | null;
  chains: Chain[];
}

// Sequence interface to group related structures
export interface RfamSequence {
  sequenceId: string;
  structures: Structure[];
  status: "already_present" | "new_sequence" | "new_structure";
  metadata: {
    length: number;
    hasBasePairs: boolean;
    hasPseudoknots: boolean;
  };
}

// RfamFamily interface for family-level data
export interface RfamFamily {
  familyId: string;
  familyName: string;
  status: "complete" | "incomplete" | "curate" | "skipped" | "no_matches";
  sequences: RfamSequence[];
  metadata: {
    totalStructures: number;
    totalSequences: number;
  };
}

// Type aliases for common status values
export type FamilyStatus = RfamFamily["status"];
export type SequenceStatus = RfamSequence["status"];
export type ChainStatus = Chain["status"];

// Utility type for filtering options
export interface RfamFilters {
  hasBasePairs: boolean;
  hasPseudoknots: boolean;
  familyStatus: "all" | FamilyStatus;
  sequenceStatus: "all" | SequenceStatus;
  chainStatus: "all" | ChainStatus;
}

// Type for summary statistics
export interface RfamSummary {
  totalFamilies: number;
  totalSequences: number;
  totalStructures: number;
  // New structures are a subset of new sequences
  familiesWithNewStructures: number;
  familiesWithNewSequences: number; // Includes both new_structure and new_sequence
  familiesWithHighResolution: number;
  statusBreakdown: {
    complete: number;
    incomplete: number;
    curate: number;
    skipped: number;
    no_matches: number;
  };
  sequenceBreakdown: {
    already_present: number;
    new_sequence: number; // Includes both new_structure and new_sequence
    new_structure: number; // Subset of new_sequence
  };
}
