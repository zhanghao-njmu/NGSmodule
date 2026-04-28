/**
 * TanStack Query hooks for the Files domain.
 */
import { useMutation, useQuery, useQueryClient } from '@tanstack/react-query'

import { queryKeys } from '@/lib/queryClient'
import fileService from '@/services/file.service'

export function useFileList(params?: any) {
  return useQuery({
    queryKey: queryKeys.files.list(params),
    queryFn: () => fileService.getAll(params),
  })
}

export function useFilesBySample(sampleId: string | undefined, params?: any) {
  return useQuery({
    queryKey: [...queryKeys.files.all, 'by-sample', sampleId, params],
    queryFn: () => fileService.getFilesBySample(sampleId as string, params),
    enabled: !!sampleId,
  })
}

export function useFilesByProject(projectId: string | undefined, params?: any) {
  return useQuery({
    queryKey: [...queryKeys.files.all, 'by-project', projectId, params],
    queryFn: () => fileService.getFilesByProject(projectId as string, params),
    enabled: !!projectId,
  })
}

export function useFile(id: string | undefined) {
  return useQuery({
    queryKey: queryKeys.files.detail(id ?? ''),
    queryFn: () => fileService.getById(id as string),
    enabled: !!id,
  })
}

export function useUploadFile() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: (params: { sample_id: string; file: File; onProgress?: (p: number) => void }) =>
      fileService.uploadFile(params),
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: queryKeys.files.all })
      queryClient.invalidateQueries({ queryKey: queryKeys.stats.storage })
    },
  })
}

export function useBatchUploadFiles() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: (params: { sampleId: string; files: File[]; onProgress?: (p: number) => void }) =>
      fileService.batchUpload(params.sampleId, params.files, params.onProgress),
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: queryKeys.files.all })
      queryClient.invalidateQueries({ queryKey: queryKeys.stats.storage })
    },
  })
}

export function useDeleteFile() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: (id: string) => fileService.delete(id),
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: queryKeys.files.all })
      queryClient.invalidateQueries({ queryKey: queryKeys.stats.storage })
    },
  })
}

export function useVerifyFileChecksum() {
  return useMutation({
    mutationFn: (id: string) => fileService.verifyChecksum(id),
  })
}
