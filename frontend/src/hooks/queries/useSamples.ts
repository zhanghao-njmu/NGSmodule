/**
 * TanStack Query hooks for the Samples domain.
 *
 * The sample service is built on the CRUD factory + a few domain helpers
 * (CSV import/export, batch create). All mutations invalidate the
 * project-scoped sample list so the UI stays in sync after writes.
 */
import { useMutation, useQuery, useQueryClient } from '@tanstack/react-query'

import { queryKeys } from '@/lib/queryClient'
import sampleService from '@/services/sample.service'

export function useSampleList(projectId: string | undefined, params?: any) {
  return useQuery({
    queryKey: queryKeys.samples.list(projectId ?? '', params),
    queryFn: () =>
      projectId
        ? sampleService.getSamplesByProject(projectId, params)
        : sampleService.getAll(params),
    enabled: !!projectId,
  })
}

export function useSample(id: string | undefined) {
  return useQuery({
    queryKey: queryKeys.samples.detail(id ?? ''),
    queryFn: () => sampleService.getById(id as string),
    enabled: !!id,
  })
}

// ----- mutations -----

export function useCreateSample() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: (data: any) => sampleService.create(data),
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: queryKeys.samples.all })
      queryClient.invalidateQueries({ queryKey: queryKeys.stats.samples })
    },
  })
}

export function useUpdateSample() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: ({ id, data }: { id: string; data: any }) => sampleService.update(id, data),
    onSuccess: (_, { id }) => {
      queryClient.invalidateQueries({ queryKey: queryKeys.samples.detail(id) })
      queryClient.invalidateQueries({ queryKey: queryKeys.samples.all })
    },
  })
}

export function useDeleteSample() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: (id: string) => sampleService.delete(id),
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: queryKeys.samples.all })
      queryClient.invalidateQueries({ queryKey: queryKeys.stats.samples })
    },
  })
}

export function useImportSamplesFromCSV() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: ({ projectId, file }: { projectId: string; file: File }) =>
      sampleService.importFromCSV(projectId, file),
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: queryKeys.samples.all })
      queryClient.invalidateQueries({ queryKey: queryKeys.stats.samples })
    },
  })
}

export function useBatchCreateSamples() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: (data: any) => sampleService.batchCreateSamples(data),
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: queryKeys.samples.all })
      queryClient.invalidateQueries({ queryKey: queryKeys.stats.samples })
    },
  })
}
