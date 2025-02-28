// dataService.ts - Updated to correctly handle chainId format

import {
  FamilyReport,
  CandidateReport,
  ChainStatus,
  FamilyStatus,
  SequenceStatus,
  RfamFamily,
  RfamSequence,
  Structure,
  Chain,
  PdbChainId,
} from "../types";

// Mock data for development environment
import mockData from "../data/rfam-data-complete.json";

export class RfamService {
  private static instance: RfamService;
  private apiUrl =
    "https://ftp.ebi.ac.uk/pub/databases/Rfam/.preview/3d/report.txt"; // Default production API endpoint
  private isDev = process.env.NODE_ENV === "development";

  private constructor() {}

  public static getInstance(): RfamService {
    if (!RfamService.instance) {
      RfamService.instance = new RfamService();
    }
    return RfamService.instance;
  }

  /**
   * Fetches raw data from the API and transforms it into the dashboard format
   */
  public async getData(): Promise<RfamFamily[]> {
    try {
      const rawData = await this.fetchFromAPI();
      const transformedData = this.transformData(rawData);
      return transformedData;
    } catch (error) {
      console.error("[RfamService] Error in getData():", error);
      throw error;
    }
  }

  /**
   * Handles fetching from either local mock data (development) or API (production)
   */
  private async fetchFromAPI(): Promise<FamilyReport[]> {
    // Use mock data in development mode
    if (this.isDev) {
      console.log("[RfamService] Using mock data in development mode");
      try {
        if (mockData) {
          console.log(
            "[RfamService] Mock data length:",
            Array.isArray(mockData) ? mockData.length : "Not an array",
          );
        }
        return Promise.resolve(mockData as FamilyReport[]);
      } catch (error) {
        console.error("[RfamService] Error loading mock data:", error);
        throw new Error(`Failed to load mock data: ${error}`);
      }
    }

    // Fetch from API in production mode
    try {
      const response = await fetch(this.apiUrl);

      if (!response.ok) {
        throw new Error(`Failed to fetch data: ${response.statusText}`);
      }

      const data = await response.json();
      return data;
    } catch (error) {
      throw error;
    }
  }

  /**
   * Transforms the raw API data into the dashboard-friendly format
   */
  private transformData(rawData: FamilyReport[]): RfamFamily[] {
    return rawData.map((family) => this.transformFamily(family));
  }

  /**
   * Extracts a standardized chainId string from various possible formats
   * @param candidateChainId The chainId from a candidate (can be string or object)
   * @returns A chainId in the format "pdbId_chainId"
   */
  private extractChainId(candidateChainId: any): {
    pdbId: string;
    chainId: string;
    fullChainId: string;
  } {
    let pdbId = "unknown";
    let chainId = "unknown";
    let fullChainId = "";

    // If it's already a string in pdb_chain format
    if (typeof candidateChainId === "string") {
      if (candidateChainId.includes("_")) {
        const parts = candidateChainId.split("_");
        pdbId = parts[0];
        chainId = parts.slice(1).join("_"); // In case chain id contains underscores
        fullChainId = candidateChainId;
      } else {
        // If it doesn't have the expected format, use as is
        pdbId = "unknown";
        chainId = candidateChainId;
        fullChainId = `${pdbId}_${chainId}`;
      }
    }
    // If it's a complex object with nested pdbId
    else if (
      typeof candidateChainId === "object" &&
      candidateChainId !== null
    ) {
      if (candidateChainId.pdbId && candidateChainId.chainId) {
        // Handle nested pdbId object
        if (
          typeof candidateChainId.pdbId === "object" &&
          candidateChainId.pdbId.pdbId
        ) {
          pdbId = candidateChainId.pdbId.pdbId;
        } else if (typeof candidateChainId.pdbId === "string") {
          pdbId = candidateChainId.pdbId;
        }

        chainId = candidateChainId.chainId;
        fullChainId = `${pdbId}_${chainId}`;
      }
    }

    return { pdbId, chainId, fullChainId };
  }

  /**
   * Transforms a FamilyReport into an RfamFamily
   */
  private transformFamily(family: FamilyReport): RfamFamily {
    // Group candidates by sequence ID
    const sequenceGroups = new Map<string, CandidateReport[]>();

    // First pass: group candidates by sequence ID
    family.candidates.forEach((candidate) => {
      if (candidate.sequence && candidate.chain) {
        // Extract chainId information
        const chainIdInfo = this.extractChainId(candidate.chainId);

        const seqId =
          candidate.sequence.sequenceId || `unknown-${chainIdInfo.fullChainId}`;
        if (!sequenceGroups.has(seqId)) {
          sequenceGroups.set(seqId, []);
        }
        sequenceGroups.get(seqId)!.push(candidate);
      }
    });

    // Transform sequence groups into RfamSequence objects
    const sequences: RfamSequence[] = [];

    sequenceGroups.forEach((candidates, sequenceId) => {
      // Group structures by PDB ID
      const structureGroups = new Map<string, CandidateReport[]>();

      candidates.forEach((candidate) => {
        // Extract chainId information
        const chainIdInfo = this.extractChainId(candidate.chainId);

        if (!structureGroups.has(chainIdInfo.pdbId)) {
          structureGroups.set(chainIdInfo.pdbId, []);
        }
        structureGroups.get(chainIdInfo.pdbId)!.push(candidate);
      });

      // Convert to Structure objects
      const structures: Structure[] = Array.from(structureGroups.entries()).map(
        ([pdbId, candidates]) => {
          const chains: Chain[] = candidates.map((candidate) => {
            // Extract chainId information
            const chainIdInfo = this.extractChainId(candidate.chainId);

            // Use the chain object to get properties (with appropriate fallbacks)
            return {
              chainId: chainIdInfo.chainId,
              status: candidate.chain
                ? this.mapChainStatus(candidate.chain.status)
                : "error",
              hasBasePairs: candidate.chain
                ? !!candidate.chain.hasBasePairs
                : false,
              hasPseudoknots: candidate.chain
                ? !!candidate.chain.hasPseudoknots
                : false,
              pdbId: chainIdInfo.pdbId,
            };
          });

          return {
            pdbId,
            method: "", // This info might not be available in the new format
            resolution: -1, // This info might not be available in the new format
            chains,
          };
        },
      );

      // Determine sequence status based on chain statuses
      const hasNewStructure = candidates.some(
        (c) => c.chain!.status === ChainStatus.NEW_CHAIN,
      );
      const status = hasNewStructure
        ? "new_structure"
        : candidates[0].sequence!.status === SequenceStatus.NEW_SEQUENCE
          ? "new_sequence"
          : "committed";

      // Create sequence object
      const sequence: RfamSequence = {
        sequenceId: sequenceId,
        structures,
        status,
        metadata: {
          length: candidates[0].sequence!.sequenceLength,
          hasBasePairs: candidates.some((c) => c.chain!.hasBasePairs),
          hasPseudoknots: candidates.some((c) => c.chain!.hasPseudoknots),
        },
      };

      sequences.push(sequence);
    });

    // Additional logic to handle cases where there are candidates without sequences
    // (e.g., errors, skipped)
    family.candidates
      .filter((c) => !c.sequence)
      .forEach((candidate) => {
        if (!candidate.chain) return; // Skip completely empty candidates

        // Extract chainId information
        const chainIdInfo = this.extractChainId(candidate.chainId);

        const chain: Chain = {
          chainId: chainIdInfo.chainId,
          status: this.mapChainStatus(candidate.chain.status),
          hasBasePairs: candidate.chain.hasBasePairs,
          hasPseudoknots: candidate.chain.hasPseudoknots,
          pdbId: chainIdInfo.pdbId,
        };

        // Create a "placeholder" sequence for error/skipped chains
        const sequence: RfamSequence = {
          sequenceId: `error-${chainIdInfo.fullChainId}`,
          structures: [
            {
              pdbId: chainIdInfo.pdbId,
              method: "",
              resolution: -1,
              chains: [chain],
            },
          ],
          status: "error",
          metadata: {
            length: 0,
            hasBasePairs: candidate.chain.hasBasePairs,
            hasPseudoknots: candidate.chain.hasPseudoknots,
          },
        };

        sequences.push(sequence);
      });

    // Determine best resolution across all structures if available
    const bestResolution =
      sequences
        .flatMap((seq) => seq.structures)
        .filter((s) => s.resolution > 0)
        .reduce((best, s) => Math.min(best, s.resolution), Infinity) || -1;

    // Create the family object
    return {
      familyId: family.familyId,
      familyName: family.familyName,
      status: this.mapFamilyStatus(family.status),
      sequences,
      metadata: {
        totalStructures: sequences.reduce(
          (sum, seq) => sum + seq.structures.length,
          0,
        ),
        totalSequences: sequences.length,
        bestResolution,
      },
    };
  }

  /**
   * Maps ChainStatus enum to display-friendly string
   */
  private mapChainStatus(status: ChainStatus): string {
    const statusMap: Record<ChainStatus, string> = {
      [ChainStatus.COMMITTED_CHAIN]: "committed",
      [ChainStatus.NEW_CHAIN]: "new_structure",
      [ChainStatus.ERROR]: "error",
      [ChainStatus.SKIPPED_CHAIN]: "skipped",
    };
    return statusMap[status];
  }

  /**
   * Maps FamilyStatus enum to display-friendly string
   */
  private mapFamilyStatus(status: FamilyStatus): string {
    const statusMap: Record<FamilyStatus, string> = {
      [FamilyStatus.COMPLETE]: "complete",
      [FamilyStatus.CURATE]: "curatable",
      [FamilyStatus.INCOMPLETE]: "updated",
      [FamilyStatus.NO_VALID]: "no_matches",
      [FamilyStatus.SKIPPED]: "skipped",
    };
    return statusMap[status];
  }
}
