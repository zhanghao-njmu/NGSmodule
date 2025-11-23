/**
 * File Store - Global state management for files
 */
import { create } from 'zustand'
import { fileService } from '../services/file.service'
import type { FileItem } from '../types/file'
import { message } from 'antd'

interface FileStore {
  // State
  files: FileItem[]
  currentFile: FileItem | null
  loading: boolean
  error: string | null
  uploadProgress: number

  // Actions
  fetchFiles: (params?: { sample_id?: string; project_id?: string; file_type?: string }) => Promise<void>
  fetchFileById: (id: string) => Promise<void>
  uploadFile: (sampleId: string, file: File) => Promise<FileItem | null>
  downloadFile: (id: string, filename: string) => Promise<void>
  deleteFile: (id: string) => Promise<void>
  setCurrentFile: (file: FileItem | null) => void
  clearError: () => void
}

export const useFileStore = create<FileStore>((set) => ({
  // Initial state
  files: [],
  currentFile: null,
  loading: false,
  error: null,
  uploadProgress: 0,

  // Fetch all files
  fetchFiles: async (params) => {
    set({ loading: true, error: null })
    try {
      const response = await fileService.getAll(params)
      set({ files: response.items, loading: false })
    } catch (error: any) {
      set({ error: error.message, loading: false })
      message.error(`Failed to fetch files: ${error.message}`)
    }
  },

  // Fetch file by ID
  fetchFileById: async (id) => {
    set({ loading: true, error: null })
    try {
      const file = await fileService.getById(id)
      set({ currentFile: file, loading: false })
    } catch (error: any) {
      set({ error: error.message, loading: false })
      message.error(`Failed to fetch file: ${error.message}`)
    }
  },

  // Upload file
  uploadFile: async (sampleId, file) => {
    set({ loading: true, error: null, uploadProgress: 0 })
    try {
      const uploadedFile = await fileService.uploadFile({
        sample_id: sampleId,
        file,
        onProgress: (percent) => {
          set({ uploadProgress: percent })
        },
      })

      set((state) => ({
        files: [uploadedFile, ...state.files],
        loading: false,
        uploadProgress: 0,
      }))

      message.success('File uploaded successfully')
      return uploadedFile
    } catch (error: any) {
      set({ error: error.message, loading: false, uploadProgress: 0 })
      message.error(`Failed to upload file: ${error.message}`)
      return null
    }
  },

  // Download file
  downloadFile: async (id, filename) => {
    set({ loading: true, error: null })
    try {
      await fileService.downloadFile(id, filename)
      set({ loading: false })
      message.success('File download started')
    } catch (error: any) {
      set({ error: error.message, loading: false })
      message.error(`Failed to download file: ${error.message}`)
    }
  },

  // Delete file
  deleteFile: async (id) => {
    set({ loading: true, error: null })
    try {
      await fileService.delete(id)
      set((state) => ({
        files: state.files.filter((f) => f.id !== id),
        currentFile: state.currentFile?.id === id ? null : state.currentFile,
        loading: false,
      }))
      message.success('File deleted successfully')
    } catch (error: any) {
      set({ error: error.message, loading: false })
      message.error(`Failed to delete file: ${error.message}`)
    }
  },

  // Set current file
  setCurrentFile: (file) => {
    set({ currentFile: file })
  },

  // Clear error
  clearError: () => {
    set({ error: null })
  },
}))

export default useFileStore
