// Updated types.ts with correct chainId format

// Status enums matching the Python implementation
export enum SequenceStatus {
  NEW_SEQUENCE = "new_sequence",
  COMMITTED_SEQUENCE = "committed_sequence",
}

export enum BasepairingStatus {
  NEW_STRUCTURE = "new_pairing",
  COMMITTED_STRUCTURE = "committed_pairing",
}

export enum MatchStatus {
  AUTOMATED_MATCH = "automated_match",
  MANUALLY_FORCED = "manually_forced",
}

export enum PkStatus {
  COMMITTED_PK = "committed_pk",
  NEW_PK = "new_pk",
  NO_PK = "no_pk",
}

export enum ChainStatus {
  COMMITTED_CHAIN = "committed_chain",
  ERROR = "error",
  NEW_CHAIN = "new_chain",
  SKIPPED_CHAIN = "skipped_chain",
}

export enum FamilyStatus {
  COMPLETE = "complete",
  CURATE = "curate",
  INCOMPLETE = "incomplete",
  NO_VALID = "no_matches",
  SKIPPED = "skipped",
}

// The main type interfaces that match the Python classes
export interface ChainReport {
  chainId: string; // Format: pdb_id_chain_id (e.g. "4u3n_4")
  status: ChainStatus;
  basepairingStatus: BasepairingStatus | null;
  dotBracket: string;
  hasBasePairs: boolean;
  hasPseudoknots: boolean;
  newBasepairCount: number | null;
  removedBasepairCount: number | null;
}

export interface MatchReport {
  sequenceStartPosition: number;
  sequenceStopPosition: number;
  bitScore: number;
  eValue: number;
  cmStartPosition: number;
  cmEndPosition: number;
  hexColor: string;
  isSignificant: boolean;
  matchStatus: MatchStatus;
}

export interface SequenceReport {
  sequenceId: string | null;
  sequenceLength: number;
  sequenceMd5: string;
  status: SequenceStatus;
}

export interface PdbChainId {
  pdbId: {
    pdbId: string;
  };
  chainId: string;
}

export interface CandidateReport {
  chainId: PdbChainId | string; // Can be either a string ("pdb_id_chain_id") or a PdbChainId object
  matchReport: MatchReport;
  chain: ChainReport | null;
  sequence: SequenceReport | null;
}

export interface AlignmentReport {
  numColumnsBefore: number;
  numColumnsAfter: number;
}

export interface FamilyReport {
  familyId: string;
  familyName: string;
  rnaType: string;
  status: FamilyStatus;
  alignmentInfo: AlignmentReport | null;
  candidates: CandidateReport[];
}

// Dashboard data-specific interfaces that convert from API format to UI-friendly format
export interface Chain {
  chainId: string; // For UI display, this is just the chain identifier (e.g., "A", "4")
  status: string; // Using string for more flexibility in display
  hasBasePairs: boolean;
  hasPseudoknots: boolean;
  pdbId: string; // PDB identifier (e.g., "4u3n")
}

export interface Structure {
  pdbId: string;
  method: string;
  resolution: number;
  chains: Chain[];
}

export interface RfamSequence {
  sequenceId: string;
  structures: Structure[];
  status: string;
  metadata: {
    length: number;
    hasBasePairs: boolean;
    hasPseudoknots: boolean;
  };
}

export interface RfamFamily {
  familyId: string;
  familyName: string;
  status: string;
  sequences: RfamSequence[];
  metadata: {
    totalStructures: number;
    totalSequences: number;
    bestResolution: number;
  };
}

// Utility type for filtering options
export interface RfamFilters {
  hasBasePairs: boolean;
  hasPseudoknots: boolean;
  familyStatus: string;
  sequenceStatus: string;
  chainStatus: string;
}

// Type for summary statistics
export interface RfamSummary {
  totalFamilies: number;
  totalSequences: number;
  totalStructures: number;
  familiesWithNewStructures: number;
  familiesWithNewSequences: number;
  familiesWithHighResolution: number;
  statusBreakdown: Record<string, number>;
  sequenceBreakdown: Record<string, number>;
}
