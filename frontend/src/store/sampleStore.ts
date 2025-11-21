/**
 * Sample Store - Global state management for samples
 */
import { create } from 'zustand'
import { sampleService } from '../services/sample.service'
import type { Sample, SampleCreate, SampleUpdate, SampleBatchCreate } from '../types/sample'
import { message } from 'antd'

interface SampleStore {
  // State
  samples: Sample[]
  currentSample: Sample | null
  loading: boolean
  error: string | null

  // Actions
  fetchSamples: (params?: { project_id?: string }) => Promise<void>
  fetchSampleById: (id: string) => Promise<void>
  createSample: (data: SampleCreate) => Promise<Sample | null>
  batchCreateSamples: (data: SampleBatchCreate) => Promise<void>
  importFromCSV: (projectId: string, file: File) => Promise<void>
  updateSample: (id: string, data: SampleUpdate) => Promise<void>
  deleteSample: (id: string) => Promise<void>
  setCurrentSample: (sample: Sample | null) => void
  clearError: () => void
}

export const useSampleStore = create<SampleStore>((set) => ({
  // Initial state
  samples: [],
  currentSample: null,
  loading: false,
  error: null,

  // Fetch all samples
  fetchSamples: async (params) => {
    set({ loading: true, error: null })
    try {
      const response = await sampleService.getSamples(params)
      set({ samples: response.items, loading: false })
    } catch (error: any) {
      set({ error: error.message, loading: false })
      message.error(`Failed to fetch samples: ${error.message}`)
    }
  },

  // Fetch sample by ID
  fetchSampleById: async (id) => {
    set({ loading: true, error: null })
    try {
      const sample = await sampleService.getSample(id)
      set({ currentSample: sample, loading: false })
    } catch (error: any) {
      set({ error: error.message, loading: false })
      message.error(`Failed to fetch sample: ${error.message}`)
    }
  },

  // Create new sample
  createSample: async (data) => {
    set({ loading: true, error: null })
    try {
      const sample = await sampleService.createSample(data)
      set((state) => ({
        samples: [sample, ...state.samples],
        loading: false,
      }))
      message.success('Sample created successfully')
      return sample
    } catch (error: any) {
      set({ error: error.message, loading: false })
      message.error(`Failed to create sample: ${error.message}`)
      return null
    }
  },

  // Batch create samples
  batchCreateSamples: async (data) => {
    set({ loading: true, error: null })
    try {
      const response = await sampleService.batchCreateSamples(data)
      message.success(`${response.count} samples created successfully`)
      // Refresh samples list
      await sampleService.getSamples({ project_id: data.project_id }).then((res) => {
        set({ samples: res.items, loading: false })
      })
    } catch (error: any) {
      set({ error: error.message, loading: false })
      message.error(`Failed to create samples: ${error.message}`)
    }
  },

  // Import from CSV
  importFromCSV: async (projectId, file) => {
    set({ loading: true, error: null })
    try {
      const response = await sampleService.importFromCSV(projectId, file)
      message.success(`${response.count} samples imported successfully`)
      // Refresh samples list
      await sampleService.getSamples({ project_id: projectId }).then((res) => {
        set({ samples: res.items, loading: false })
      })
    } catch (error: any) {
      set({ error: error.message, loading: false })
      message.error(`Failed to import samples: ${error.message}`)
    }
  },

  // Update sample
  updateSample: async (id, data) => {
    set({ loading: true, error: null })
    try {
      const updatedSample = await sampleService.updateSample(id, data)
      set((state) => ({
        samples: state.samples.map((s) => (s.id === id ? updatedSample : s)),
        currentSample:
          state.currentSample?.id === id ? updatedSample : state.currentSample,
        loading: false,
      }))
      message.success('Sample updated successfully')
    } catch (error: any) {
      set({ error: error.message, loading: false })
      message.error(`Failed to update sample: ${error.message}`)
    }
  },

  // Delete sample
  deleteSample: async (id) => {
    set({ loading: true, error: null })
    try {
      await sampleService.deleteSample(id)
      set((state) => ({
        samples: state.samples.filter((s) => s.id !== id),
        currentSample: state.currentSample?.id === id ? null : state.currentSample,
        loading: false,
      }))
      message.success('Sample deleted successfully')
    } catch (error: any) {
      set({ error: error.message, loading: false })
      message.error(`Failed to delete sample: ${error.message}`)
    }
  },

  // Set current sample
  setCurrentSample: (sample) => {
    set({ currentSample: sample })
  },

  // Clear error
  clearError: () => {
    set({ error: null })
  },
}))

export default useSampleStore
