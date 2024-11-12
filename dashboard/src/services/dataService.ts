import type { RfamFamily, RfamSequence, Structure, Chain } from "../types";
import localData from "../data/rfam-data.json";

function isValidChain(chain: unknown): chain is Chain {
  return (
    typeof chain === "object" &&
    chain !== null &&
    "chainId" in chain &&
    "status" in chain &&
    "hasBasePairs" in chain &&
    "hasPseudoknots" in chain &&
    typeof (chain as Chain).chainId === "string" &&
    typeof (chain as Chain).hasBasePairs === "boolean" &&
    typeof (chain as Chain).hasPseudoknots === "boolean" &&
    [
      "already_present",
      "new_sequence",
      "new_structure",
      "error",
      "skipped",
    ].includes((chain as Chain).status)
  );
}

export function isValidStructure(structure: unknown): structure is Structure {
  return (
    typeof structure === "object" &&
    structure !== null &&
    "pdbId" in structure &&
    "method" in structure &&
    "resolution" in structure &&
    "chains" in structure &&
    typeof (structure as Structure).pdbId === "string" &&
    ((structure as Structure).method === null ||
      typeof (structure as Structure).method === "string") &&
    ((structure as Structure).resolution === null ||
      typeof (structure as Structure).resolution === "number") &&
    Array.isArray((structure as Structure).chains) &&
    (structure as Structure).chains.every(isValidChain)
  );
}

function isValidSequence(sequence: unknown): sequence is RfamSequence {
  return (
    typeof sequence === "object" &&
    sequence !== null &&
    "sequenceId" in sequence &&
    "structures" in sequence &&
    "status" in sequence &&
    "metadata" in sequence &&
    Array.isArray((sequence as RfamSequence).structures) &&
    (sequence as RfamSequence).structures.every(isValidStructure) &&
    ["already_present", "new_sequence", "new_structure"].includes(
      (sequence as RfamSequence).status,
    )
  );
}

export function isValidRfamFamily(family: unknown): family is RfamFamily {
  return (
    typeof family === "object" &&
    family !== null &&
    "familyId" in family &&
    "familyName" in family &&
    "status" in family &&
    "sequences" in family &&
    "metadata" in family &&
    typeof (family as RfamFamily).familyId === "string" &&
    typeof (family as RfamFamily).familyName === "string" &&
    ["complete", "incomplete", "curate", "skipped", "no_matches"].includes(
      (family as RfamFamily).status,
    ) &&
    Array.isArray((family as RfamFamily).sequences) &&
    (family as RfamFamily).sequences.every(isValidSequence)
  );
}

export class RfamServiceError extends Error {
  constructor(
    message: string,
    public readonly statusCode?: number,
    public readonly originalError?: unknown,
  ) {
    super(message);
    this.name = "RfamServiceError";
  }
}

export class RfamService {
  private static instance: RfamService;
  private readonly apiUrl: string;
  private readonly dataSource: string;

  private constructor() {
    this.apiUrl = import.meta.env.VITE_API_URL;
    this.dataSource = import.meta.env.VITE_DATA_SOURCE;
  }

  public static getInstance(): RfamService {
    if (!RfamService.instance) {
      RfamService.instance = new RfamService();
    }
    return RfamService.instance;
  }

  private groupStructuresBySequence(family: any): RfamSequence[] {
    // Create a map to group structures by sequenceId
    const sequenceMap = new Map<
      string,
      {
        structures: Structure[];
        hasBasePairs: boolean;
        hasPseudoknots: boolean;
        allChainsNewStructure: boolean;
        hasNewStructure: boolean;
      }
    >();

    // First pass: collect all structures for each unique sequenceId
    family.structures.forEach((structure: Structure) => {
      structure.chains.forEach((chain) => {
        if (chain.sequenceId) {
          const existingSequence = sequenceMap.get(chain.sequenceId);
          if (!existingSequence) {
            sequenceMap.set(chain.sequenceId, {
              structures: [structure],
              hasBasePairs: chain.hasBasePairs,
              hasPseudoknots: chain.hasPseudoknots,
              allChainsNewStructure: chain.status === "new_structure",
              hasNewStructure: chain.status === "new_structure",
            });
          } else {
            if (
              !existingSequence.structures.find(
                (s) => s.pdbId === structure.pdbId,
              )
            ) {
              existingSequence.structures.push(structure);
            }
            existingSequence.hasBasePairs =
              existingSequence.hasBasePairs || chain.hasBasePairs;
            existingSequence.hasPseudoknots =
              existingSequence.hasPseudoknots || chain.hasPseudoknots;
            existingSequence.allChainsNewStructure =
              existingSequence.allChainsNewStructure &&
              chain.status === "new_structure";
            existingSequence.hasNewStructure =
              existingSequence.hasNewStructure ||
              chain.status === "new_structure";
          }
        }
      });
    });

    // Convert map to array of RfamSequences
    return Array.from(sequenceMap.entries()).map(([sequenceId, data]) => ({
      sequenceId,
      structures: data.structures,
      // If any chain is new_structure, mark as new_structure
      // Otherwise, it's a new_sequence (since it's in our data)
      status: data.hasNewStructure ? "new_structure" : "new_sequence",
      metadata: {
        length: 0, // This would need to be populated with actual sequence length if available
        hasBasePairs: data.hasBasePairs,
        hasPseudoknots: data.hasPseudoknots,
      },
    }));
  }

  private calculateFamilyStats(family: any) {
    const stats = {
      totalStructures: family.structures.length,
      totalNewStructures: 0,
      totalNewSequenceStructures: 0,
      hasNewStructures: false,
      hasNewSequences: false,
      hasBasePairs: false,
      hasPseudoknots: false,
    };

    // Count structures based on their chain statuses
    family.structures.forEach((structure: Structure) => {
      const hasNewStructure = structure.chains.some(
        (chain) => chain.status === "new_structure",
      );
      const hasNewSequence = structure.chains.some(
        (chain) =>
          chain.status === "new_structure" || chain.status === "new_sequence",
      );

      if (hasNewStructure) {
        stats.totalNewStructures++;
        stats.hasNewStructures = true;
        stats.hasNewSequences = true; // new_structure implies new_sequence
      } else if (hasNewSequence) {
        stats.totalNewSequenceStructures++;
        stats.hasNewSequences = true;
      }

      // Check for structural features
      structure.chains.forEach((chain) => {
        if (chain.hasBasePairs) stats.hasBasePairs = true;
        if (chain.hasPseudoknots) stats.hasPseudoknots = true;
      });
    });

    return stats;
  }

  private transformLegacyData(families: any[]): RfamFamily[] {
    return families.map((family) => {
      const sequences = this.groupStructuresBySequence(family);
      const stats = this.calculateFamilyStats(family);

      return {
        familyId: family.familyId,
        familyName: family.familyName,
        status: family.status,
        sequences,
        metadata: {
          totalStructures: stats.totalStructures,
          totalSequences: sequences.length,
          totalNewStructures: stats.totalNewStructures,
          totalNewSequenceStructures: stats.totalNewSequenceStructures,
          hasNewStructures: stats.hasNewStructures,
          hasNewSequences: stats.hasNewSequences,
          hasBasePairs: stats.hasBasePairs,
          hasPseudoknots: stats.hasPseudoknots,
        },
      };
    });
  }

  private sortSequencesByStatus(families: RfamFamily[]): RfamFamily[] {
    return families.map((family) => ({
      ...family,
      sequences: [...family.sequences].sort((a, b) => {
        const statusPriority = {
          new_structure: 0,
          new_sequence: 1,
          already_present: 2,
        };
        return statusPriority[a.status] - statusPriority[b.status];
      }),
    }));
  }

  private async fetchFromApi(): Promise<RfamFamily[]> {
    try {
      const response = await fetch(this.apiUrl);
      if (!response.ok) {
        throw new RfamServiceError(
          `HTTP error! status: ${response.status}`,
          response.status,
        );
      }
      const data = await response.json();
      if (!Array.isArray(data)) {
        throw new RfamServiceError("Invalid data structure received from API");
      }

      const transformedData = this.transformLegacyData(data);

      if (!transformedData.every(isValidRfamFamily)) {
        throw new RfamServiceError(
          "Invalid data structure after transformation",
        );
      }

      return this.sortSequencesByStatus(transformedData);
    } catch (error) {
      if (error instanceof RfamServiceError) {
        throw error;
      }
      throw new RfamServiceError("Failed to fetch Rfam data", undefined, error);
    }
  }

  private getLocalData(): RfamFamily[] {
    const transformedData = this.transformLegacyData(localData);

    if (!transformedData.every(isValidRfamFamily)) {
      throw new RfamServiceError("Invalid local data structure");
    }

    return this.sortSequencesByStatus(transformedData);
  }

  public async getData(): Promise<RfamFamily[]> {
    if (this.dataSource === "local") {
      return this.getLocalData();
    }
    return this.fetchFromApi();
  }
}
